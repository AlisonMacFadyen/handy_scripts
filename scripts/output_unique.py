# Basic script to compare values in two columns and output unique values.

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help=".csv containing columns you want to compare.")
parser.add_argument("-c1", "--column1", help="Header for column you want to find unique values in.")
parser.add_argument("-c2", "--column2", help="Header for column you want to compare.")
parser.add_argument("-o", "--output", help=".csv to save unique values to.")
args = parser.parse_args()


# Read the CSV file
df = pd.read_csv(args.input)

# Find values in column X that are not in column Y
unique_to_x = df[~df[args.column1].isin(df[args.column2])][args.column1].unique()

# Print the unique values
# print("Values unique to column X:")
# for value in unique_to_x:
#     print(value)

# Save the unique values to a new CSV file
pd.DataFrame(unique_to_x, columns=['Unique_Values']).to_csv(args.output, index=False)