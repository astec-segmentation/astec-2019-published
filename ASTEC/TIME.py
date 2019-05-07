import os
import os.path

def append_lines_to_file(filename, lines_to_append):
    if os.path.exists(filename)==True:
        with open(filename, 'a') as file:
            for l in lines_to_append:
                file.write(str(l)+"\n")
    else:
        with open(filename, 'w') as file:
            for l in lines_to_append:
                file.write(str(l)+"\n")