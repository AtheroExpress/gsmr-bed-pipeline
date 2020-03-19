import configparser

config = configparser.ConfigParser()
config.read('pipeline.config')
section = config['DEFAULT']

for key in section:
    globals()[key] = section[key]
