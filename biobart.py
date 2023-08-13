import os
import re
import pandas as pd
import argparse
import torch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def cli():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(prog='BioBART', description='BioBART v1.0 -- Generate short description for InterProScan functional annotations.')

    # Add an argument for the input file path
    parser.add_argument('-i', '--input', type=str, help='Path to inpterpro tsv file', action='store', required=True, metavar='')
    parser.add_argument('-f', '--fasta', type=str, help='Path to fasta file', action='store', required=True, metavar='')
    parser.add_argument('-l', '--long_text', help='Optional, Set the length of the generated text. Default is one liner output, if this option is included, it would generate 2 to 3 descriptive sentences.', 
                        action='store_true', default=False)
    parser.add_argument('-o', '--output-file', type=str, help='Optional, Output file name. Please do not add any file extension, as the output would be in fasta format', action='store', default='ProAnno', metavar='')
    parser.add_argument('-v', '--version', help='Optional, Display the version of BioBART.', action='version', version='%(prog)s 1.0')

    # Parse the command line arguments
    args = parser.parse_args()
    action = ProAnno(**vars(args))
    action.output()

class ProAnno:
    def __init__(self, input='store', fasta = 'store', long_text = False, output_file = 'store'):
        self.input = input
        self.fasta = fasta
        self.longtext = long_text
        self.output_file = output_file

        if not os.path.splitext(input)[1].lower() == '.tsv':
            raise TypeError('Input file for -i, interpro must be a tsv file')
        
        if not os.path.splitext(fasta)[1].lower() == '.fasta':
            raise TypeError('Input file for -f must be a fasta file')


    def preprocess(self):
        '''Processing of the tsv and fasta file into data list with the respective name, annotations and sequence.'''
        input = pd.read_csv(self.input, sep='\t')
        records = list(SeqIO.parse(self.fasta, 'fasta'))
        df = input.drop(input.columns[1:11], axis=1)
        df.columns = range(len(df.columns))
        df = df.reset_index(drop=True)
        df = df.loc[(df.iloc[:, 1] != '-')]
        df = df.drop_duplicates()

        ## ---Modified from chatgpt generated code---
        data_list = df.groupby(0).agg({2: list}).reset_index().values.tolist()

        for data in data_list:
            for record in records:
                if record.id == data[0]:
                    data.append(str(record.seq))
                    break

        return data_list

        ## ---Modified from chatgpt generated code---

    def word_count(self, sentence):
        '''Calcualtion of word count within an input text.'''
        return len(sentence.split())
    
    def output(self):
        '''Generation of funnctional annotation passage and summarisation for output towards fasta file.'''
        from transformers import BioGptTokenizer, BioGptForCausalLM, set_seed
        from transformers import AutoTokenizer, BartForConditionalGeneration
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # BioGPT
        model_1 = BioGptForCausalLM.from_pretrained("microsoft/BioGPT")
        tokenizer_1 = BioGptTokenizer.from_pretrained("microsoft/BioGPT")
        model_1.to(device)

        # BART Summarize        
        model_2 = BartForConditionalGeneration.from_pretrained("facebook/bart-large-cnn")
        tokenizer_2 = AutoTokenizer.from_pretrained("facebook/bart-large-cnn")
        model_2.to(device)

        # with open(f'{self.output_file}_BART.txt', "w") as f:
        #         f.write(f"")
        with open(f'{self.output_file}.fasta', "w") as f:
                f.write(f"")

        node_no = 0
        nodes = self.preprocess()
        for node in nodes:
            node_id = node[0]
            functions = node[1]
            functions_str = ", ".join(functions)
            sequence = node[2]
            bart_response_first_sentence = f"{node_id} is annotated with {functions_str}. "
            bart_response_str = ''
            for function in functions:
            # BioGPT Generation
                sentence = f"{function} is"
                inputs = tokenizer_1(sentence, return_tensors="pt").to(device)
                set_seed(42)
                with torch.no_grad():
                    beam_output = model_1.generate(**inputs,
                                            min_length=100,
                                            max_length=1024,
                                            num_beams=5,
                                            early_stopping=True
                                            ).to(device)
                    bio_output = tokenizer_1.decode(beam_output[0], skip_special_tokens=True)

                # For BART Summarisation about each annotation
                inputs = tokenizer_2([bio_output], max_length=1024, truncation=True, return_tensors="pt").to(device)
                summary_ids = model_2.generate(inputs["input_ids"], num_beams=2, min_length=30, max_length=200).to(device)
                bart_response = tokenizer_2.batch_decode(summary_ids, skip_special_tokens=True, clean_up_tokenization_spaces=False)[0] 

                bart_response_str += bart_response + ' '
            
            # BART final summarisation of all annotation passage if word cound is more than 50 words    
            if self.word_count(bart_response_str) >= 50:    
                inputs = tokenizer_2([bart_response_str], max_length=1024, truncation=True, return_tensors="pt").to(device)
                summary_ids = model_2.generate(inputs["input_ids"], num_beams=2, min_length=30, max_length=100).to(device)
                bart_second_response = tokenizer_2.batch_decode(summary_ids, skip_special_tokens=True, clean_up_tokenization_spaces=False)[0]
            else:
                bart_second_response = bart_response_str

            # To split and only use the first sentence by default and more if -l is chosen
            if self.longtext == False:
                 final_bart_response = re.split(r'(?<=[.!?])\s', bart_second_response)[0]
            else:
                 final_bart_response = bart_second_response

            w_output = SeqRecord(Seq(sequence), id = node_id, description= str(bart_response_first_sentence + final_bart_response))
            with open(f'{self.output_file}.fasta', "a") as output_handle:
                    SeqIO.write(w_output, output_handle,'fasta')
            node_no +=1
            
        output_handle.close()
        torch.cuda.empty_cache
        print('BioBART Run Completed.')

if __name__ == '__main__':
     cli()
