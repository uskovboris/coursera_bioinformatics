def read_file(file_path):
    """ read and return specified file content """
    with open(file_path) as f:
        text = f.read()
        f.close()
    return text


def read_lines(file_path):
    """ read file as a list of lines """
    with open(file_path) as f:
        return [l.strip() for l in f.readlines()]


def write_file(file_path, text):
    """ rewrite file content """
    with open(file_path, "w") as f:
        f.write(text)
