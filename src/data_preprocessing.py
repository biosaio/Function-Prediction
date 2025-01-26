import pandas as pd
import numpy as np
import h5py
from Bio import SeqIO
import tensorflow as tf
import csv

def load_fasta(config):
  sequences = {}
  with open(config, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
      sequences[record.id] = str(record.seq)
  return sequences

def load_go_embeddings(config):
  go_embeddings = {}
  with open(config['paths']['go_embeddings_csv'], mode='r') as f:
      reader = csv.reader(f)
      header = next(reader)

      for row in reader:
          go_term = row[0]
          embedding = row[1]
          embedding_list = [float(x) for x in embedding.strip('[]').split(',')]
          go_embeddings[go_term] = embedding_list

  return go_embeddings

def load_train_data(config):
  # Train set
  train_set = pd.read_csv(config['paths']['train_set'], sep='\t')
  train_set.rename(columns={'Protein_ID': 'protein_id', 'GO_term': 'go_term'}, inplace=True)
  train_set.sort_values(by=['protein_id'], inplace=True)

  # Train sequences
  train_sequences = load_fasta(config['paths']['train_fasta'])

  # Protein embeddings
  with h5py.File(config['paths']['train_embeddings'], 'r') as f:

    # Train IDs
    train_ids = np.array(list(f.keys()))
    train_embeddings = np.array([f[train_id][:] for train_id in train_ids])

  # InterPro domains
  train_domains = pd.read_csv(config['paths']['train_protein2ipr'], sep='\t', header=None, names=config['parameters']['DOMAINS_COLS'])

  return train_set, train_sequences, train_ids, train_embeddings, train_domains

def load_test_data(config):
  # Test sequences
  test_sequences = load_fasta(config['paths']['test_fasta'])

  # Protein embeddings
  with h5py.File(config['paths']['test_embeddings'], 'r') as f:

    # Test IDs
    test_ids = np.array(list(f.keys()))
    test_embeddings = np.array([f[test_id][:] for test_id in test_ids])

  # InterPro domains
  test_domains = pd.read_csv(config['paths']['test_protein2ipr'], sep='\t', header=None, names=config['parameters']['DOMAINS_COLS'])

  # BLAST results
  test_blast = pd.read_csv(config['paths']['test_blast'], sep='\t')

  return test_sequences, test_ids, test_embeddings, test_domains, test_blast

def add_sequence_embedding(sequence_embeddings, df):
  """
  sequence_embeddings: ndarray of embedded protein sequences
  df: dataframe with protein IDs
  """
  sequence_embeddings_dict = {protein_id: sequence_embedding for protein_id, sequence_embedding in zip(df['protein_id'], sequence_embeddings)}
  df['sequence_embedding'] = df['protein_id'].map(sequence_embeddings_dict)

  return df

def add_sequence_length(sequences_dict, df):
  """
  sequences_dict: dictionary of sequences
  df: dataframe with protein IDs
  """
  sequences_lengths = {protein_id: len(sequence) for protein_id, sequence in sequences_dict.items()}
  df['sequence_length'] = df['protein_id'].map(sequences_lengths)

  return df

def add_go_embedding(go_embeddings, df):
  """
  go_embeddings: dictionary of GO terms and their embeddings
  df: dataframe with GO terms
  """
  df['go_embedding'] = df['go_term'].map(go_embeddings)

  return df

def add_domain(domains_df, df):
  """
  domains_df: dataframe with InterPro domains
  df: dataframe with protein IDs
  """
  domains_dict = domains_df.groupby('protein_id')['ipr_id'].unique().to_dict()
  df['domains'] = df['protein_id'].map(domains_dict)

  return df

def explode_by_domain(df):
  df = df.explode('domains')
  df.rename(columns={'domains': 'ipr_id'}, inplace=True)

  return df

def split_by_aspect(df):
  """
  df: a dataframe with GO aspect columns
  """
  df_molecular_function = df[df['aspect'] == 'molecular_function']
  df_biological_process = df[df['aspect'] == 'biological_process']
  df_cellular_component = df[df['aspect'] == 'cellular_component']

  return df_molecular_function, df_biological_process, df_cellular_component
