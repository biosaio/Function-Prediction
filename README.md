# Function-Prediction
Predicting functional GO annotation given Prot5 embeddings and BLAST results.

Steps for training: 
- blast the training data against itself to obtain proteins that are similar to each other
- preprocess the result, to build a table where we count how many times a GO term appears in the list of similar proteins for each protein
- apply TF-IDF transformation to this table
- train Multi Layer Perceptrons with ProtT5 embeddings concatenated with our TF-IDF table. We will have 3 models, one for each GO aspect

To evaluate:
- blast test set against the train set
- apply TF-IDF transformation with coefficients calculated in training
- use the models to predict functions
- (optional, if you have ground-truth): use CAFAEVAL to evaluate the results

A first step, produce the BLAST results using the standalone blast program for linux command line. 
Create blast dataset with train sequences:

`makeblastdb -in train.fasta -parse_seqids -dbtype prot -out blastdb `

Break into chunks:

`seqkit split -s 5000 train.fasta -O train_chunks`

Then use the scripts process_chunks.sh and process_chunks_rm.sh to obtain the blast results and clean them, to be ready to be used in the python pipeline.

Use preprocessBlast notebook to preprocess the BLAST chunks, obtaining count tables for GO terms associated with the proteins significant blast hits. 
Then use IF_IDF_preprocessing_final to finalize the preprocessing, creating the tsv files ready to be fed to the training procedure. Here we produce the TF-IDF transformed blast count tables and concatenate them with 
the protein embeddings
Finally use the train notebook to train the models and evaluate the performance on the test set using cafaeval. 
