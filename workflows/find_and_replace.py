import os


def replace_text_in_files(directory, old_text, new_text):
    # Walk through all files in the specified directory
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            # Read the current contents of the file
            try:
                with open(filepath, 'r', encoding='utf-8') as file:
                    filedata = file.read()

                # Replace 'wdl_common' with 'wdl_common'
                new_filedata = filedata.replace(old_text, new_text)

                # Write the changes back to the file only if a change occurred
                if new_filedata != filedata:
                    with open(filepath, 'w', encoding='utf-8') as file:
                        file.write(new_filedata)
            except Exception as e:
                print(f"Failed to process {filepath}: {e}")

# Replace 'wdl_common' with 'wdl_common' in all .wdl files in the current directory
directory = os.getcwd()
replace_text_in_files(directory, 'wdl_common', 'wdl_common')
