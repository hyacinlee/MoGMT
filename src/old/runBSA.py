import argparse

def main(args=None):
    parser = argparse.ArgumentParser(description="function1 的帮助信息")
    parser.add_argument("-s", "--source", help="源参数")
    parser.add_argument("-v", "--value", help="值参数")

    # 如果 args 为 None，则自动解析 sys.argv
    parsed_args = parser.parse_args(args)

    print(f"function1 被执行: -s {parsed_args.source} -v {parsed_args.value}")

if __name__ == "__main__":
    main() 
