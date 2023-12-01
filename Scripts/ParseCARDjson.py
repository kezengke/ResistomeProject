import pandas as pd

json_file_path = "Andermann_pipeline_01312023/REFs/CARD/card.json"

df = pd.read_json(json_file_path)

num_columns = df.shape[1]
# Delete the last three columns
df = df.drop(df.columns[num_columns - 3:], axis=1)

# Get the number of rows in the DataFrame
num_rows = df.shape[0]

# Delete the last three rows
df = df.drop(df.index[num_rows - 2:])
df_transposed = df.transpose()

sorted_df = df_transposed.sort_values(by='model_id')

sorted_df.to_csv("CARD_geneLength.csv", index=False)
