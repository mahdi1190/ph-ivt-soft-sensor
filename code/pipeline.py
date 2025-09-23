"""Example CLI pipeline to fetch/process data and reproduce results."""
from pathlib import Path
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['fetch','process','figures'], required=True)
    p.add_argument('--item', default=None, help='Dataset or figure ID')
    args = p.parse_args()
    Path('figures').mkdir(exist_ok=True)
    print(f'Pretend to run stage={args.stage}, item={args.item}')

if __name__ == "__main__":
    main()