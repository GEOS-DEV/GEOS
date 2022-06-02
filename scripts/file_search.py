import sys, getopt
import os
import glob
import subprocess


def html_head():
    txt = '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">'
    txt += "<head><style>body{font-size:90%;font-family: Arial, Helvetica, sans-serif;padding:5em;} \
    #script {font-family:monospace; padding-left:5em;} emph{color:red;}</style></head><body>"
    return txt


def html_tail():
    txt = "</body>"
    return txt


def collect_files(root_folder, ext):
    """
    Crawls in all directories from a list of folders
    and returns all .ext files in any subfolders
    """
    abs_input_file_list = []
    for root, dirs, files in os.walk(root_folder):
        for file in files:
            if file.endswith(ext):
                abs_input_file_list.append(os.path.join(root, file))

    print("Found " + str(len(abs_input_file_list)) + " " + ext + " files")
    return abs_input_file_list


def tokenize(text, token, highlight=""):
    text = text.replace("--\n--\n", "</br>")
    text = text.replace("\n", "</br>")
    if highlight:
        text = text.replace(highlight, "<emph>{}</emph>".format(highlight))

    return "<{0:}>{1:}</{0:}>".format(token, text)


def parse_and_search(allfiles, todo_token, n_lines_before=5, n_lines_after=5):
    print("Looking for {} in {} files...".format(todo_token, len(allfiles)))
    cnt = 0
    log_output = "output.html"
    with open(log_output, "w") as f:
        f.write(html_head())
        f.write(tokenize("Instances of {}".format(todo_token), "h1"))
        f.write(tokenize("Number of lines before : {}".format(n_lines_before), "p"))
        f.write(tokenize("Number of lines after  : {}".format(n_lines_after), "p"))
        for file_name in allfiles:
            p = subprocess.Popen(
                [
                    "egrep",
                    "-i",
                    "-n",
                    "-A{}".format(n_lines_after),
                    "-B{}".format(n_lines_before),
                    todo_token,
                    file_name,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            outlog, errlog = p.communicate()
            output_log = outlog.decode("UTF-8")
            if len(output_log) > 0:
                cnt += 1
                f.write(tokenize(file_name, "h4"))
                f.write(tokenize(output_log, "div id=script", todo_token))
        summary = "</br>Found {} files containing {} (case insensitive) ".format(
            cnt, todo_token.upper()
        )
        f.write(tokenize(summary, "strong"))
        f.write(html_tail())

    print(summary)
    print("Results are written to file: output.html")


def main(argv):
    folder_to_search = os.getcwd()
    token = "TODO"
    try:
        opts, args = getopt.getopt(argv, "hf:t:", ["folder_to_search=", "token="])
    except getopt.GetoptError:
        print("test.py -f <folder_to_search> -t <token>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":
            print("file_search.py -f <folder_to_search> -t <token>")
            sys.exit()
        elif opt in ("-f", "--folder_to_search"):
            folder_to_search = os.path.abspath(arg)
        elif opt in ("-t", "--token"):
            token = arg

    if not (os.path.exists(folder_to_search)):
        print("*** Folder does not exist : {}".format(folder_to_search))
        sys.exit()
    print("Folder to search is : {}".format(folder_to_search))
    print("Token  to search is : {}".format(token))

    hppfiles = collect_files(folder_to_search, "hpp")
    cppfiles = collect_files(folder_to_search, "cpp")
    all_files = hppfiles + cppfiles
    parse_and_search(all_files, token)


if __name__ == "__main__":
    main(sys.argv[1:])
