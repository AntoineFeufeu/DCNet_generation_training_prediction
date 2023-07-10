###################
## This file is used to generate indexes for train -test-split.
##Notice: Just execute 1 time, because each time it will give diffent result randomly
##################
import logging
import argparse
from numpy import savetxt
from Feed_Data import DataGenerator
def read_args():
    parser = argparse.ArgumentParser()                   
    parser.add_argument("--total",
                        default=842,
                        type=int,
                        help="total number of samples")
    parser.add_argument("--deb",
                        default=0,
                        type=int,
                        help="total number of samples")
    parser.add_argument("--ratio",
                        default=0.8,
                        type=float,
                        help="ratio of split")
    args = parser.parse_args()
    return args
def main(args):
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    data_generator = DataGenerator()
    train_idx, valid_idx, test_idx = data_generator.generate_split_indexes(args.deb,args.total,args.ratio)
    savetxt('./train_idx.csv', train_idx, delimiter=',')
    savetxt('./valid_idx.csv', valid_idx, delimiter=',')
    savetxt('./test_idx.csv', test_idx, delimiter=',')
    return

if __name__ == '__main__':
  args = read_args()
  main(args)
