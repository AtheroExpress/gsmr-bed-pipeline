import configparser

config = configparser.ConfigParser(strict=False)
config.read('pipeline.config')

def load_section(section='DEFAULT'):
    section = config[section]

    for key in section:
        globals()[key] = section[key]
