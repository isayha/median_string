# CPSC 450 - Bioinformatics
# Isayha Raposo - 230133508
# Assignment 1 - Median String Problem

from sys import argv, exit
from os import path
from math import inf
from itertools import product

bases = ['A', 'T', 'C', 'G']

# Returns all possible k-length combinations of the four bases (A, T, C, G)
def get_possible_k_mers(k):
    k_mers = list(product(bases, repeat = k))
    return k_mers

# Returns the Hamming distance between the provided k-mer and the provided (sub)string (see minimum_distance())
def distance(k_mer, string):
    # Exits for cases in which the length of the provided (sub)string is not equal to the length of the provided k-mer
    if len(k_mer) != len(string):
        print("ERROR: len(k_mer) != len(substring)")
        print("DEBUGGING REQUIRED")
        exit(1)

    distance = 0
    for index in range(0, len(k_mer)):
        if k_mer[index] != string[index]:
            distance += 1
    return distance

# Returns the MINIMUM Hamming distance between the provided k-mer and the provided string
# This is achieved by dividing the string into substrings of the same length as the k-mer and passing the two to minimum_distance()
def minimum_distance(k_mer, string):
    minimum_distance = inf
    substring_start = 0
    substring_end = len(k_mer)
    while substring_end <= len(string):
        substring = string[substring_start:substring_end]
        current_distance = distance(k_mer, substring)
        if current_distance < minimum_distance:
            minimum_distance = current_distance
        substring_start += 1
        substring_end += 1
    return minimum_distance

# The algorithm (as defined via pseudocode in Chapter 2)
def median_string(k, DNA_seqs):
    k_mers = get_possible_k_mers(k)
    k_mer_min_dist_sums = []
    for k_mer in k_mers:
        sum_of_min_dists = 0
        for DNA_seq in DNA_seqs:
            sum_of_min_dists += minimum_distance(k_mer, DNA_seq)
        k_mer_min_dist_sums.append((k_mer, sum_of_min_dists))
    median_string = min(k_mer_min_dist_sums, key = lambda x : x[1])[0]
    return median_string 

# Returns the (int) k and (list) DNA sequences found within the specified/provided data file
def process_data_file(data_file_name):
    data_file = open(data_file_name, 'r')
    k = int(data_file.readline())

    if k < 1 or k > 10:
        print("ERROR: Invalid k value specified. See README.md.")
        exit(0)

    DNA_seqs = []
    data_file_line = data_file.readline()

    for DNA_seq in data_file_line.split(' '): # Handles space delimiting
        # The assignment details state that DNA sequences should be newline-separated BUT the test data provided is SPACE-separated
        # This function can handle both means of delimiting
        for base in DNA_seq.upper():
            if base not in bases:
                print("ERROR: Invalid DNA sequence/nucleotide/base specified.")
                exit(0)

        while DNA_seq:
            if len(DNA_seq) < k:
                print("ERROR: Specified DNA sequence has length less than specified k value.")
                exit(0)
            DNA_seqs.append(DNA_seq.upper())
            DNA_seq = data_file.readline()

    return k, DNA_seqs

# Prints the provided string to the console and writes it to the provided output_file
def print_and_write(output_file, string):
    print(string)
    output_file.write(string + '\n')

def main():
    arg_count = len(argv)
    if arg_count != 2:
        print("ERROR: No command line argument provided. See README.md.")
        exit(1)
    data_file_name = argv[1]
    if not path.isfile(data_file_name):
        print("ERROR: Specified data file " + data_file_name + " not found. See README.md.")
        exit(1)
    
    k, DNA_seqs = process_data_file(data_file_name)

    print("Specified data file found and processed:", "\nk =", k, "\nDNA sequences:", DNA_seqs, '\n')

    output_file_name = data_file_name.split(".txt")[0] + "_output.txt"
    output_file = open(output_file_name, "w")
    
    solution = median_string(k, DNA_seqs)
    formatted_solution = ''
    for base in solution:
        formatted_solution += base
    print_and_write(output_file, "Solution: " + formatted_solution)
    print_and_write(output_file, "(Note that more than one solution may exist: In such cases, only one solution is provided by this program as per the assignment details)")
    print("\nThis solution has been written to " + output_file_name)

if __name__ == "__main__":
    main()