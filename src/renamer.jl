using FileIO
# using StringCase

# Define the directory path
main_directory_path  = "./motifs/totalenergy/Japan"

# Define the word to capitalize
word_to_capitalize = "japan"

# Get a list of all files in the directory
files = readdir(directory_path)

# Function to capitalize the first letter of a specific word in a string
function capitalize_word_in_string(s::AbstractString, word::AbstractString)
    parts = split(s, "_")
    for (i, part) in enumerate(parts)
        if part == word
            parts[i] = uppercasefirst(part)
        end
    end
    return join(parts, "_")
end

# # Iterate through each file and rename it
# for file in files
#     # Check if it's a file (not a directory)
#     if isfile(joinpath(directory_path, file))
#         # Capitalize the specified word in the file name
#         new_name = capitalize_word_in_string(file, word_to_capitalize)
        
#         # Construct the full paths of the old and new file names
#         old_path = joinpath(directory_path, file)
#         new_path = joinpath(directory_path, new_name)
        
#         # Rename the file
#         mv(old_path, new_path)
        
#         println("Renamed '$file' to '$new_name'")
#     end
# end

# Recursive function to iterate through files in subdirectories
function process_files_in_directory(directory_path::AbstractString)
    # Get a list of all items (files and subdirectories) in the directory
    items = readdir(directory_path)
    
    for item in items
        item_path = joinpath(directory_path, item)
        
        if isfile(item_path)
            # If it's a file, capitalize the specified word in the filename
            new_name = capitalize_word_in_string(item, word_to_capitalize)
            
            # Construct the full paths of the old and new file names
            old_path = item_path
            new_path = joinpath(directory_path, new_name)
            
            # Rename the file
            mv(old_path, new_path)
            
            println("Renamed '$item' to '$new_name'")
        elseif isdir(item_path)
            # If it's a subdirectory, recursively process its files
            process_files_in_directory(item_path)
        end
    end
end

# Start processing files in the main directory (including subdirectories)
process_files_in_directory(main_directory_path)