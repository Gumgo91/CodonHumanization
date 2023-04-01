<h1 align="center">
  CodonHumanization: A Novel Approach to Codon Optimization Based on Machine Learning and Genetic Algorithm 
  <br/>
</h1>

## Table of Contents

* [Abstract](#Abstract)
* [Installation](#Installation)
* [Method](#Method)
* [Contact](#Contact)

## Abstract
Codon optimization is a molecular biological technique that modulates codon usage and modifies DNA or RNA sequence to enhance protein expression. However, for human, current codon optimization methods have certain limitations, such as potential alterations in protein structure and increased immunogenicity. In this study, we present a new approach to codon optimization utilizing machine learning and genetic algorithm, termed "codon humanization". The utilization of natural language processing (NLP) can enable the detection of latent patterns in human DNA sequences, which can subsequently be employed for the humanization of sequences through directed evolution using genetic algorithms. This approach has the potential to improve safety and efficacy in gene delivery to humans, such as gene therapy or RNA vaccines.

Codon humanization can be performed using an already trained classification model. In addition, all materials for learning and reproduction of genetic algorithms are provided.

## Installation


## Method
For the purpose of implementing codon humanization, a dataset was compiled and a classification model was developed. The binary classification model can differentiate between a random sequence and a human sequence when presented with a DNA sequence. By employing directed evolution via genetic algorithms on a DNA sequence, the corresponding amino acid sequence remains unchanged, but the DNA sequence is transformed to become more human-like. The figure below shows the sequential processes in which a binary classification model is trained, and a given sequence is humanized using a genetic algorithm.
![Figure 1](https://user-images.githubusercontent.com/65825773/229262827-ee810488-bdb4-44d0-8559-72354946a5ac.png)
