# main.py
import argparse
import dataset as data
import functions as hp
import importlib
import numpy as np
from sympy import pprint
importlib.reload(hp)
importlib.reload(data)


def main():
    parser = argparse.ArgumentParser(
        description="Computations for the paper: knotoid colorings by biquandles."
    )

    parser.add_argument(
        "action",
        choices=["count","bvb_matrix","validate","all"],
        help="Action to perform."
    )

    parser.add_argument(
        "--code",
        required=False,
        help="Knotoid code label (e.g. 3.1.10)."
    )

    parser.add_argument(
        "--biquandle",
        required=False,
        type=int,
        help="Biquandle index (integer: 0, 1, 2, ...)."
    )

    parser.add_argument(
        "--bvb",
        required=False,
        type=int,
        help="Biquandle virtual bracket index (integer: 0, 1, 2, ...)."
    )

    parser.add_argument(
        "--list",
        action="store_true",
        help="List available virtual knotoid labeled peer codes and biquandle indices, then exit."
    )

    args = parser.parse_args()

    # List mode
    if args.list:
        print("Available virtual knotoid labeled peer codes :")
        print("  " + ", ".join(sorted(c.label for c in data.parsed_codesV)))

        if args.action == "count":
            print("\nAvailable biquandle indices:")
            print("  " + ", ".join(str(i) for i in range(len(data.bq_list))))

        elif args.action == "bvb_matrix":
            print("\nAvailable biquandle virtual bracket indices:")
            print("  " + ", ".join(str(i) for i in range(len(data.bvb))))

    if args.action == "count":
        if args.code is None or args.biquandle is None:
            raise SystemExit("Usage: python main.py count --code <LABEL> --biquandle <INDEX>")

        if not (0 <= args.biquandle < len(data.bq_list)):
            raise SystemExit(
                f"Biquandle index out of range: {args.biquandle}. "
                f"Valid range: 0–{len(data.bq_list)-1}"
            )

        code_label = args.code
        bq = args.biquandle
        for code in data.parsed_codesV:
            if code.label==code_label:
                coloring_list, coloring_list_filt = hp.color(code, data.bq_list[bq].p ,data.bq_list[bq].under,data.bq_list[bq].over)
                print("The number of colorings of", code.label," with biquandle B",bq, "is :", len(coloring_list))
                if len(coloring_list)>0:
                    print("The colorings of", code.label," with biquandle X",bq, "are :")
                    for i in coloring_list_filt:
                        print(i)
                n = data.bq_list[bq].p
                counting_matrix = np.zeros((n, n), dtype=int)
                for coloring in coloring_list_filt:
                    counting_matrix[coloring[0],coloring[len(coloring)-1]]+=1
                print("X",bq,"-counting matrix of ",code.label, "is :")
                print(counting_matrix)
    elif args.action == "bvb_matrix":
        if args.code is None or args.bvb is None:
            raise SystemExit(
                "Usage: python main.py bvb_matrix --code <LABEL> --bvb <INDEX>"
            )
        code_label = args.code
        bvb = args.bvb
        matrix = hp.bvb_m(code_label,bvb)
        pprint(matrix)
    elif args.action == "validate":
        if args.bvb is None and args.biquandle is None:
            raise SystemExit(
                "Usage: python main.py validate --biquandle <INDEX> or python main.py validate --bvb <INDEX>"
            )
        if args.biquandle is None:
            i = args.bvb
            under = data.bvb[i].under
            over = data.bvb[i].over
            p = data.bvb[i].bvbp
            A = data.bvb[i].A
            B = data.bvb[i].B
            V = data.bvb[i].V
            C = data.bvb[i].C
            D = data.bvb[i].D
            U = data.bvb[i].U
            all_good, delta, omega = hp.is_biquandle_bracket(under,over,p,A,B,V,C,D,U)
            if all_good == True:
                print("This is a valid Z_",len(under),"virtual bracket over the ring Z_",data.bvb[i].bvbp)
                print("Delta =",delta)
                print("Omega =",omega)
        elif args.bvb is None:
            if hp.is_biquandle(data.bq_list[args.biquandle].under,data.bq_list[args.biquandle].over):
                print("These under/over operation matrices")
                pprint(data.bq_list[args.biquandle].under)
                pprint(data.bq_list[args.biquandle].over)
                print("defines a biquandle structure over Z_",data.bq_list[args.biquandle].p)
    elif args.action == "all":
        if args.bvb is None and args.biquandle is None:
            raise SystemExit(
                "Usage: python main.py validate --biquandle <INDEX> or python main.py validate --bvb <INDEX>"
            )
        if args.biquandle is None:
            for code in data.parsed_codesV:
                code_label = code.label
                bvb = args.bvb
                matrix = hp.bvb_m(code_label,bvb)
                print("Biquandle virtual bracket matrix of ",code.label,"w.r.t biquandle virtual bracket",bvb,"is :")
                pprint(matrix)
        if args.bvb is None:
            bq = args.biquandle
            for code in data.parsed_codesV:
                coloring_list, coloring_list_filt = hp.color(code, data.bq_list[bq].p ,data.bq_list[bq].under,data.bq_list[bq].over)
                n = data.bq_list[bq].p
                counting_matrix = np.zeros((n, n), dtype=int)
                for coloring in coloring_list_filt:
                    counting_matrix[coloring[0],coloring[len(coloring)-1]]+=1
                print("Biquandle counting matrix of ",code.label,"with biquandle",bq,"is :")
                print(counting_matrix)

if __name__ == "__main__":
    main()
