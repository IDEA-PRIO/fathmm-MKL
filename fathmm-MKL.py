#!/usr/bin/env python

import argparse
import subprocess
from multiprocessing.pool import ThreadPool

'''
fathmm-MKL.py: Predict the Functional Consequences of Single Nucleotide Variants (SNVs)
'''

# fetch argument(s)
parser = argparse.ArgumentParser(description='Predict the Functional Consequences of Single Nucleotide Variants (SNVs)',
                                 add_help=False)
parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
parser.add_argument("-threads", type=int, default=1)

group = parser.add_argument_group("Required")
group.add_argument('input_file', metavar='<F1>', type=argparse.FileType("r"), help='the mutation data to process')
group.add_argument('output_file', metavar='<F2>', type=argparse.FileType("w"), help='where predictions are written')
group.add_argument('db', metavar='<db>', type=argparse.FileType("r"),
                   help='precomputed database of fathmm-MKL predictions')
args = parser.parse_args()


def build_queries():
    queries = []
    errors = []
    for line in args.input_file:
        if not line.strip() or line.startswith("#"):
            continue
        line_parts = line.strip().upper().split(",")
        try:
            assert line_parts.__len__() == 4  # required data present in query
            int(line_parts[1])  # is position numeric
            assert line_parts[2] in ["A", "C", "G", "T"]  # expected base
            assert line_parts[3] in ["A", "C", "G", "T"]  # expected base
        except:
            errors.append("\t".join(
                ['', '', '', '', '', '', '', '', "Error: Unexpected Format '" + ",".join(line_parts) + "'"]) + "\n")
            continue
        query_string = f"{line_parts[0]}:{int(line_parts[1]) + 1}-{int(line_parts[1]) + 1}"
        query = (query_string, line_parts)
        queries.append(query)
    return queries, errors


def run_query(query):
    query_string, query_parts = query
    query_result = ['', '', '', '', "No Prediction Found"]
    command = f"tabix {args.db.name} {query_string}"
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    data, err = proc.communicate()

    if err:
        query_result[-1] = "Error: 'tabix' command";
    elif data:
        for record in data.decode().split("\n"):
            if not record:
                continue
            record = record.strip().split("\t")

            if not record[0] == query_parts[0]:
                query_result[-1] = "Error: Unexpected Chromosome";
                break
            if not record[1] == query_parts[1]:
                query_result[-1] = "Error: Unexpected Position";
                break
            if not record[3] == query_parts[2]:
                query_result[-1] = "Warning: Inconsistent Base (Expecting '" + record[3] + "')";
                break
            if record[4] == query_parts[3]:
                query_result = record[5:] + ['']
                break
    return (query_result, query_parts)


def run_queries():
    results = []
    queries, errors = build_queries()
    results += errors
    with ThreadPool(processes=args.threads) as pool:
        jobs = [pool.apply_async(run_query, (query,)) for query in queries]
        pool.close()
        pool.join()
        query_results = [job.get() for job in jobs]
    return results + query_results


def write_results(results):
    file_header = "\t".join(
        ["# Chromosome", "Position", "Ref. Base", "Mutant Base", "Non-Coding Score", "Non-Coding Groups",
         "Coding Score", "Coding Groups", "Warning"])
    otherlines = ["\t".join(query + result) for (result, query) in results]
    for line in [file_header] + otherlines:
        args.output_file.write(line + "\n")


def main():
    results = run_queries()
    write_results(results)


if __name__ == '__main__':
    main()
