
###############
## Libraries ##
###############

import os
import pandas as pd
from Bio import SeqIO
import re
import gzip
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import argparse
import sys

####################
## Main functions ##
####################

def get_adapt_data(fasta_file):
    adapt_seqs = []
    for s in SeqIO.parse(fasta_file, "fasta"):
        adapt_seqs.append(str(s.seq))
    print(color_text(f"{len(adapt_seqs)} adapters on file {fasta_file}", "yellow"))
    return adapt_seqs

def trim_adaptors_advanced(records, adaptor_list, percent_threshold=0.2, min_len_match=5):
    """Trims adaptor sequences or parts of adaptors from SeqRecord objects
    based on their position (start, end, or middle) and length percentage.

    This is a generator function.
    The 'records' argument should be a list or iterator returning SeqRecord objects.
    The 'adaptor_list' argument should be a list of adaptor sequences (list of str)
    to be searched for in each record.
    'percent_threshold': If an adaptor (or part) is within this percentage from an end,
                         it's considered at the start/end. Otherwise, it's in the middle.
    'min_len_match': Minimum length of a match to be considered an adaptor part.
    """
    
    for record in records:
        original_seq_str = str(record.seq)
        len_original_seq = len(original_seq_str)
        
        # Define os limites de 20%
        start_threshold_pos = int(len_original_seq * percent_threshold)
        end_threshold_pos = len_original_seq - int(len_original_seq * percent_threshold)

        best_match = None  # Armazenará a melhor correspondência encontrada (match object do re)
        best_adaptor_len = 0 # Armazenará o comprimento do adaptador que gerou o best_match
        
        # Itera sobre cada adaptador na lista fornecida
        for adaptor in adaptor_list:
            # Para cada adaptador, tenta encontrar todas as ocorrências
            # re.escape() é usado para tratar o adaptador como uma string literal,
            # caso ele contenha caracteres especiais de regex (ex: . , * + ?)
            for match in re.finditer(re.escape(adaptor), original_seq_str):
                match_start = match.start()
                match_end = match.end()
                match_len = match_end - match_start

                # Considera apenas correspondências que sejam de tamanho mínimo
                if match_len < min_len_match:
                    continue

                # Lógica para priorizar o maior match
                if match_len > best_adaptor_len:
                    best_adaptor_len = match_len
                    best_match = match
        
        trimmed_seq_str = original_seq_str # Inicializa com a sequência original
        
        # Variáveis para as qualidades e anotações
        trimmed_qualities = None
        new_letter_annotations = record.letter_annotations.copy() # Copia todas as anotações
        
        if best_match: # Se encontramos algum adaptador
            match_start = best_match.start()
            match_end = best_match.end()
            
            # Condição 1: No começo (<20%)
            if match_start < start_threshold_pos:
                trimmed_seq_str = original_seq_str[match_end:]
                # Fatia as qualidades também
                if 'phred_quality' in record.letter_annotations:
                    trimmed_qualities = record.letter_annotations['phred_quality'][match_end:]
            
            # Condição 2: No final (>80%, ou seja, a menos de 20% do final)
            elif match_end > end_threshold_pos:
                trimmed_seq_str = original_seq_str[:match_start]
                # Fatia as qualidades também
                if 'phred_quality' in record.letter_annotations:
                    trimmed_qualities = record.letter_annotations['phred_quality'][:match_start]
            
            # Condição 3: No meio (acima de 20% do começo E a menos de 20% do final)
            else:
                trimmed_seq_str = original_seq_str[:match_start] + original_seq_str[match_end:]
                # Concatena as qualidades também
                if 'phred_quality' in record.letter_annotations:
                    trimmed_qualities = record.letter_annotations['phred_quality'][:match_start] + \
                                        record.letter_annotations['phred_quality'][match_end:]
        
        # Atualiza as anotações de qualidade se existirem e foram cortadas
        if trimmed_qualities is not None:
            new_letter_annotations['phred_quality'] = trimmed_qualities
        elif 'phred_quality' in new_letter_annotations and best_match:
            # Se não houve corte, mas o adaptador foi encontrado (e o match.len era pequeno demais, por exemplo)
            # ou se a sequencia ficou vazia (o que não deveria acontecer com o slicing acima)
            # ou se o adaptador foi removido e a seq ficou vazia, mas as qualidades nao
            # Neste cenário, se não houve trimmed_qualities, e havia qualidade, o problema é que a
            # sequência cortada pode ter 0 de comprimento.
            # Vamos garantir que as qualidades correspondam ao comprimento da nova sequência
            if len(trimmed_seq_str) == 0:
                new_letter_annotations['phred_quality'] = []
            elif len(new_letter_annotations['phred_quality']) != len(trimmed_seq_str):
                 # Isso pode acontecer se houver um bug ou edge case na logica de corte acima.
                 # Em geral, o slicing de qualidades deve corresponder ao da seq.
                 # Mas como fallback, se os comprimentos não batem, remove as qualidades.
                 # A Biopython prefere não ter qualidades a ter qualidades de tamanho incorreto.
                del new_letter_annotations['phred_quality'] # Remove para evitar erro.
                print(f"Atenção: Qualidades não correspondem ao comprimento da sequência após corte para {record.id}. Qualidades removidas.")


        # Cria um novo SeqRecord com a sequência potencialmente cortada/concatenada
        # E AGORA, COM AS ANOTAÇÕES DE LETRA (incluindo qualidades) COPIADAS/AJUSTADAS!
        new_record = record.__class__(Seq(trimmed_seq_str), 
                                      id=record.id, 
                                      description=record.description,
                                      letter_annotations=new_letter_annotations)
        # Se você usa `seq.quality` diretamente em vez de `letter_annotations['phred_quality']`,
        # o Biopython pode tentar criar isso automaticamente.
        # Garantir que a propriedade quality seja definida, se necessário.
        if hasattr(record.seq, 'quality') and trimmed_qualities is not None:
            # Esta linha garante que a propriedade .quality no objeto Seq está definida
            # e pode ser importante para como SeqIO.write FASTQ lida com isso.
            new_record.seq.quality = trimmed_qualities


        yield new_record

def running_adapter_filter(fastq_file, adapters, output_file_path, pct_thr):
    if fastq_file.endswith(".gz"):
        with gzip.open(fastq_file,"rt") as handle:
            record = SeqIO.parse(handle, "fastq")

            trimmed_seqs = trim_adaptors_advanced(record, adapters, percent_threshold=pct_thr)
            SeqIO.write(trimmed_seqs, output_file_path, "fastq")
    else:
        record = SeqIO.parse(fastq_file, "fastq")

        trimmed_seqs = trim_adaptors_advanced(record, adapters, percent_threshold=pct_thr)
        SeqIO.write(trimmed_seqs, output_file_path, "fastq")

########################
## Secondary funtions ##
########################

tty_colors = {
	'green' : '\033[0;32m%s\033[0m',
	'yellow' : '\033[0;33m%s\033[0m',
	'red' : '\033[0;31m%s\033[0m'
}

# # or, example if wanting the ones from before with background highlighting
# tty_colors = {
#     'green' : '\x1b[6;37;42m%s\x1b[0m',
#     'red' : '\x1b[0;37;41m%s\x1b[0m'
# }

#Color texto de acordo com o "Warning"

def color_text(text, color='green'):

	if sys.stdout.isatty():
		return tty_colors[color] % text
	else:
		return text

################################
### Help and argument parser ###
################################

arg_parser = argparse.ArgumentParser(description = "Script for trimming reads with bad linked adapters", 
	epilog = "Ex. usage: infernommatic.py -i R1.fastq.gz -I R2.fastq.gz -o filtered_R1.fastq -O filtered_R2.fastq")


arg_parser.add_argument("-i", "--input_r1", help = "R1 fastq file -- either gzipped or not", required=True)
arg_parser.add_argument("-I", "--input_r2", help = "R2 fastq file -- either gzipped or not", required=True)
arg_parser.add_argument("-adapt", "--adapters_file", help = "Fasta file with all adapters to be processed", default="adapters.fasta")
arg_parser.add_argument("-o", "--output_r1", help = "R1 fastq output file path", required=True)
arg_parser.add_argument("-O", "--output_r2", help = "R2 fastq output file path", required=True)
arg_parser.add_argument("-pct", "--sequence_pct", help = "percent_threshold': If an adaptor (or part) is within this percentage from an end, it's considered at the start/end of the sequence.", default=0.2, type=int)

#Se nenhum comando foi dado ao script, automaticamente é mostrado o "help"

if len(sys.argv)==1:
	arg_parser.print_help(sys.stderr)
	sys.exit(0)

##################################
### Setting starting variables ###
##################################

# getting primary script full path
path = os.path.realpath(__file__)


# getting primary script directory full path
primary_script_path = path.split("/")[:-1]
primary_script_path = "/".join(primary_script_path)

args = arg_parser.parse_args()
args_dict = vars(arg_parser.parse_args())

#INPUTS
input_r1 = os.path.abspath(args_dict["input_r1"])
input_r2 = os.path.abspath(args_dict["input_r2"])

#OUTPUTS
output_r1 = os.path.abspath(args_dict["output_r1"])
output_r2 = os.path.abspath(args_dict["output_r2"])

#ADAPTERS
adapters_file = os.path.abspath(args_dict["adapters_file"])

#OPTIONS
sequence_pct = args_dict["sequence_pct"]

#####################
## Starting script ##
#####################

## Reading adapter file ##
adapters_seqs = get_adapt_data(adapters_file)

## Running filter ##
inputs = [input_r1, input_r2]
outputs = [output_r1, output_r2]

for i,o in zip(inputs, outputs):
    print(color_text(f"Running filter on file {i}"))
    running_adapter_filter(i, adapters_seqs, o, sequence_pct)
    print(color_text(f"Output saved as {o}"))

