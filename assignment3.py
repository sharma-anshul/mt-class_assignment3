#!/usr/bin/env python
import optparse
import glob
import os
import math
from collections import Counter # multiset represented by dictionary

# This code computes the AMBER extension for the BLEU metric
# to find the best correlation between refence and machine 
# translation systems.


def createNGrams(sentence, n):
  ngrams = []
  sentence = map(lambda x: x.lower(), sentence)
  for i in xrange(len(sentence) - n + 1):
    ngram = ()
    for j in xrange(n):   
      ngram += (sentence[i+j], )
    ngrams += [ngram]
  #print ngrams
  return ngrams

def score(system, ref):
  
  r = float(sum([len(x) for x in ref]))
  c = float(sum([len(x) for x in system]))

  bp = math.exp(1 - (r/c)) if c <= r else 1.0
 
  matched_unis  = sum([sum((Counter(createNGrams(s, 1)) & Counter(createNGrams(r, 1))).values()) for (s,r) in zip(system, ref)])
  uniprec = math.log(float(matched_unis) / sum([len(createNGrams(sentence, 1)) for sentence in system]))
  unirec  = math.log(float(matched_unis) / sum([len(createNGrams(sentence, 1)) for sentence in ref]))

  matched_bis   = sum([sum((Counter(createNGrams(s, 2)) & Counter(createNGrams(r, 2))).values()) for (s,r) in zip(system, ref)])
  biprec = math.log(float(matched_bis) / sum([len(createNGrams(sentence, 2)) for sentence in system]))
  birec  = math.log(float(matched_bis) / sum([len(createNGrams(sentence, 2)) for sentence in ref]))

  matched_tris  = sum([sum((Counter(createNGrams(s, 3)) & Counter(createNGrams(r, 3))).values()) for (s,r) in zip(system, ref)])
  triprec = math.log(float(matched_tris) / sum([len(createNGrams(sentence, 3)) for sentence in system]))  
  trirec  = math.log(float(matched_tris) / sum([len(createNGrams(sentence, 3)) for sentence in ref]))

  matched_quads = sum([sum((Counter(createNGrams(s, 4)) & Counter(createNGrams(r, 4))).values()) for (s,r) in zip(system, ref)])
  quadprec = math.log(float(matched_quads) / sum([len(createNGrams(sentence, 4)) for sentence in system]))
  quadrec  = math.log(float(matched_quads) / sum([len(createNGrams(sentence, 4)) for sentence in ref]))

  matched_pents = sum([sum((Counter(createNGrams(s, 5)) & Counter(createNGrams(r, 5))).values()) for (s,r) in zip(system, ref)])
  pentprec = math.log(float(matched_pents) / sum([len(createNGrams(sentence, 5)) for sentence in system]))
  pentrec  = math.log(float(matched_pents) / sum([len(createNGrams(sentence, 5)) for sentence in ref]))

  P = (uniprec + biprec + triprec + quadprec) / 4.0
  R = unirec

  alpha = 0.9

  avgP = math.pow((uniprec * biprec * triprec * quadprec), 0.25)
  Fmean = (P * R) / ((alpha * P)  + (1.0 - alpha) * (R))
  avgF = ((uniprec * unirec) / ((alpha * uniprec)  + (1.0 - alpha) * (unirec)) + (biprec * birec) / ((alpha * biprec)  + (1.0 - alpha) * (birec)) + (triprec * trirec) / ((alpha * triprec)  + (1.0 - alpha) * (trirec)) + (quadprec * quadrec) / ((alpha * quadprec)  + (1.0 - alpha) * (quadrec))) / 4.0

  return bp * (0.3 * avgP + 0.5 * Fmean + 0.2 * avgF)
 
optparser = optparse.OptionParser()
optparser.add_option("-d", "--data-dir", dest="data", default="data/dev", help="Directory containing system outputs")
(opts,_) = optparser.parse_args()
(source, reference) = (opts.data + "/source", opts.data + "/reference")

ref = [line.split() for line in open(reference)]

# read and score each system and output ordered by score
sys_scores = [(os.path.basename(f), score([line.split() for line in open(f)], ref)) for f in glob.glob(opts.data + "/*") if f != reference and f != source]
for (sysname, score) in sorted(sys_scores, key=lambda x: -x[1]):
  print sysname
