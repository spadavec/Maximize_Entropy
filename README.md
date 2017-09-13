# Maximize Entropy 

#### What in tarnation? 

This is a simple script that looks to select a subset of molecules/entries from an input csv that attempt to maximize the number of 'on' bits at different indexes within a circular fingerprints. Right now, it's really ugly. Like, really, really ugly. Put on your sunglasses and peep at the code if you dare, but it essentially does a bruteforce pairwise comparison of a "test" FP against a set of "train" FPs, and finds the one that adds the most unique 'on' bits to the collection. Once it finds that FP/molecule, it adds it to the training set, removes it from the test set, and repeats the process until you reach a threshold of 'on' bits (number of 1s / 1024). I'm sure I could be clever and try to do a multi-sort based on index, but this is quick and dirty and meant for me primarily. 

### but-why.gif 

Because it turns out that fingerprint cosine distance isn't actually the best a-priori metric for predicting  IC50/whatever-y-value-you-care-about accuracy. It turns out that 3 conditions need to be met (according to me, of course): 1) The number of 'on' bits is > 0.5, 2) Number of training examples is >~25, and 3) y-value range of training examples covers the range you actually care about. So, this script helps you find a subset that satisfy part #1. Boring, I know, but thats life, yeah? 


