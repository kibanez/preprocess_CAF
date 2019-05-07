import argparse
import ConfigParser
import logging
from CAFPreprocess import *


def read_cfg_file(config_file):
    """
    Function that reads the config file and puts together all input parameters for the main script
    :param config_file: path to the configuration file
    :return: dictionary containing parameters (as keys) and their values (as values)
    """
    file_input = open(config_file, 'r')
    config = ConfigParser.ConfigParser()
    config.readfp(file_input)

    hash_cfg = {}
    for field in config.options('INPUT'):
        hash_cfg[field] = config.get('INPUT', field)
    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT', field)
    for field in config.options('FEATURES'):
        hash_cfg[field] = config.get('FEATURES', field)
    file_input.close()
    return hash_cfg


def main():
    formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)
    logger = logging.getLogger("preprocess")
    logger.setLevel(logging.INFO)
    logger.addHandler(console)

    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", default=None, required=True,
                        help="Path to the configuration file",
                        dest="f_cfg")
    args = parser.parse_args()

    if args.f_cfg is not None:

        if not os.path.isfile(args.f_cfg):
            raise IOError('The config file %s does not exist' % args.f_cfg)

        hash_cfg = read_cfg_file(args.f_cfg)

        preprocessing_caf = CAFPreprocess(hash_cfg)

        logger.info("Checking contaminated samples...")
        preprocessing_caf.filter_out_contamination()
        logger.info("Checking reported vs estimated gender...")
        preprocessing_caf.reported_vs_genetic_gender_check()
        logger.info("Writing CAF preprocessing log file")
        f_log = preprocessing_caf.write_log()
        logger.info("Check %s to see the summary" % f_log)


if __name__ == '__main__':
    main()
