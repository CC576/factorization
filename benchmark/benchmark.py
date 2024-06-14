from fillDataset import fill
import sys
import json
import os.path
import subprocess
import resource

# 1: trial division
# 2: fermat
# 3: quadratic sieve
# 4: gnfs


def testAlg(alg: str, dataset: dict, start=0):
    outFile = "output" + alg + ".txt"
    statFile = "statistics" + alg + ".txt"

    printStats = False

    for i in range(start, len(dataset)):
        newPoint = dict(dataset[i])
        timeout = 2                # da aggiustare
        fout = open(outFile, "a")

        try:
            args = ["../build/factoring_algorithms", alg, str(newPoint["n"])]
            if(printStats): args.append("1")
            coso = subprocess.run(args, text=True, capture_output=True, timeout=timeout, check=True)

        except subprocess.TimeoutExpired:
            status = "timeout"
            print("Timeout for algrithm " + alg + " expired at numbit = ", newPoint["numbits (of n)"])
            break

        except subprocess.CalledProcessError as e:
            status = "fail"
            print("Error " + str(e.returncode) + ": ", e.stderr)
            print("Stopping experiments for algorithm " + alg)
            break

        else:
            result = coso.stdout
            result = result[1:-2].split(',')
            p = int(result[0].strip(), 10)
            q = int(result[1].strip(), 10)

            status = (p*q == newPoint["n"])
            p=abs(p); q=abs(q)
            status = status and p!=1 and p!=newPoint["n"] and q!=1 and q!= newPoint["n"]

            if(status): status = "ok"
            else: status = "math error"

            newPoint['(pRes, qRes)'] = (p, q)

            if(printStats and (alg == '3' or alg == '4')):
                stats = {"index": newPoint["index"], "nbit":newPoint["numbits (of n)"], "n":newPoint["n"], "status": status, stats:coso.stderr}

                fstat = open(statFile, "a")
                json.dump(stats, fstat)
                fstat.close()

        finally:
            newPoint['status'] = status
            # aggiungere tempo e memoria usati
            json.dump(newPoint, fout)
            fout.write(',\n')
            fout.close()





def main():
    if not os.path.isfile("dataset.json"):
        fill()

    f = open("dataset.json", "r")
    try:
        dataset = json.load(f)
    except:
        f.close()
        fill()
        f = open("dataset.json", "r")
        dataset = json.load(f)

    #print(dataset[0])
    f.close()


    if(len(sys.argv) > 1):
        alg = sys.argv[1]       # è una stringa
        start = 0
        if(len(sys.argv) > 2):
            start = int(sys.argv[2])
        testAlg(alg, dataset, start)
    else:
        for i in range(1, 5):
            testAlg(str(i), dataset)




if __name__=='__main__':
    main()