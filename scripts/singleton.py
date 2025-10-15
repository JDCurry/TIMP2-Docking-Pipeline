import pandas as pd

df = pd.read_csv("prioritized_compounds.csv")

# count how many compounds per scaffold
counts = df.groupby("scaffold").size()

# get scaffolds with exactly 1 compound
singleton_scaffolds = counts[counts == 1].index

# filter the main table for just those
singletons = df[df["scaffold"].isin(singleton_scaffolds)]

# rank them by docking_score (lowest = best)
singletons = singletons.sort_values("docking_score")

singletons.to_csv("prioritized_singletons.csv", index=False)
