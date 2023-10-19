import os
import shutil

# Function to capitalize the first letter of the word
def capitalize_first_letter(word):
    if len(word) > 0:
        return word[0].upper() + word[1:]
    else:
        return word

# Function to recursively search for files and rename them
def rename_files_with_word(root_dir, search_word):
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if search_word in file:
                # Get the index of the search word
                index = file.index(search_word)

                # Capitalize the first letter of the search word
                capitalized_search_word = capitalize_first_letter(search_word)

                # Create the new file name
                new_name = file[:index] + capitalized_search_word + file[index + len(search_word):]

                # Create the old and new file paths
                old_path = os.path.join(root, file)
                new_path = os.path.join(root, new_name)

                # Rename the file
                shutil.move(old_path, new_path)
                print(f"Renamed: {old_path} to {new_path}")

# Set the root directory and search word
root_directory = "./motifs/totalenergy/japan"  # Replace with the actual directory path
search_word = "japan"        # Replace with the word you're looking for

# Call the function to rename files
rename_files_with_word(root_directory, search_word)
