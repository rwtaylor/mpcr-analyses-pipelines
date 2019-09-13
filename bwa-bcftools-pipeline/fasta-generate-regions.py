#!/usr/bin/env python3

import sys
import re

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is copied from https://github.com/lh3/readfq
# Fast parsing of fasta/fastq
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#fasta_in    = open("/Users/ryan/Projects/testpy/test.fasta")
#regions_out = open('/Users/ryan/Projects/testpy/out.fasta', 'w+')

nmatcher      = re.compile(r'(N{1000,})') # Regex for N-blocks
regions       = [] # list of regions

# First go through each sequence and break up scaffolds into contigs. Save contig regions in a list.
# regions like this: [(length, chrom, start, end), (length, chrom, start, end), (length, chrom, start, end)]
for name, seq, qual in readfq(sys.stdin): # Iterate through each sequence in the fasta
    seq_start = 0 # reset sequence start
    if nmatcher.search(seq): # If there are any N-blocks in the sequence
        for m in nmatcher.finditer(seq): # Iterate through each N-block match
            # The end of the sequence block is 1 character before the start of the N-block
            seq_end = m.start() - 1
            # Append the sequence region prior to the start of the N-block to the regions list
            regions.append((seq_end - seq_start, name, seq_start, seq_end))
            # The next seqeunce will start 1 character after the N-block
            seq_start = m.end() + 1
        # Now that all N-blocks have been processed, add the final sequence block (if any)
        if seq_start < len(seq):
            regions.append((len(seq) - seq_start, name, seq_start, len(seq)))
    # If no N-blocks, just print the entire sequence region
    else: regions.append((len(seq) - seq_start, name, seq_start, len(seq)))
    # Before 
# Now place the regions into tasks
task_seq_limit  = 5e6 # The maximum total length of regions per task (unless region is larger than this, then task will only consist of that 1 region)
region_seq_min  = 0 # The minimum length of sequence for a region to be included
task_counter    = 0  # counter for each task
task_seq_cumsum = 0  # cumulative sum of sequences added to each task
tasks           = [] # list of tasks to be returned

for length, name, seq_start, seq_end in regions:
    if(length >= region_seq_min):
        # if adding this sequence to task will bump the task over the limit, append the task, set the new task as region
        if task_seq_cumsum > 0 and length + task_seq_cumsum > task_seq_limit:
            tasks.append(str(task_counter) + "," + str(task_seq_cumsum) + "," + task)
            task_counter += 1
            task_seq_cumsum = 0
        if task_seq_cumsum == 0:
            task = "--region " + name + ":" + str(seq_start) + "-" + str(seq_end)
            task_seq_cumsum = length
        else:
            task += " --region " + name + ":" + str(seq_start) + "-" + str(seq_end)
            task_seq_cumsum += length

for task in tasks:
    sys.stdout.write(task + "\n")
