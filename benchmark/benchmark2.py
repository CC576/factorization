#from fillDataset import fill
import sys
import json
#import os.path
import subprocess
import resource

# 1: trial division
# 2: fermat
# 3: quadratic sieve
# 4: gnfs


def testAlg(alg: str, dataset: dict, start=0, end=None):
    outFile = "output" + alg + ".txt"
    statFile = "statistics" + alg + ".txt"

    printStats = False
    timeout = 60*60*24*15       # = 15 giorni
    statsThreshold = 1800       # = 30 min

    if end is None:
        end = len(dataset)

    for i in range(start, end):
        print("current index: ", i)
        newPoint = dict(dataset[i])
        fout = open(outFile, "a")                       # se il file non esiste viene creato

        try:
            args = ["../build/factoring_algorithms", alg, str(newPoint["n"])]
            if(printStats): args.append("1")
            coso = subprocess.run(args, text=True, capture_output=True, timeout=timeout, check=True)

        except subprocess.TimeoutExpired:
            status = "timeout"
            newPoint["time"] = timeout
            print("Timeout for algrithm " + alg + " expired at numbit = ", newPoint["numbits (of n)"])
            break

        except subprocess.CalledProcessError as e:
            status = "fail"
            print("Error " + str(e.returncode) + ": ", e.stderr)
            print("Stopping experiments for algorithm " + alg)
            break

        else:
            output = coso.stdout
            #print(output)

            result = output.split('\n')[0]
            result = result[1:-1].split(',')
            #print(result)
            p = int(result[0].strip(), 10)
            q = int(result[1].strip(), 10)

            # controllare che l'output sia corretto
            status = (p*q == newPoint["n"])
            p=abs(p); q=abs(q)
            status = status and p!=1 and p!=newPoint["n"] and q!=1 and q!= newPoint["n"]

            if(status): status = "ok"
            else: status = "math error"

            newPoint['(pRes, qRes)'] = (p, q)

            # tempo e memoria usati
            resUsage = output.split('\n')[1]
            resUsage = "{" + resUsage + "}"
            #print(resUsage)
            resUsage = json.loads(resUsage)
            newPoint.update(resUsage)

            # stampare statistiche
            if(printStats and (alg == '3' or alg == '4')):
                stats = {"index": newPoint["index"], "nbit":newPoint["numbits (of n)"], "n":newPoint["n"], "status": status, "stats":coso.stderr.strip()}

                fstat = open(statFile, "a")
                json.dump(stats, fstat)
                fstat.write(',\n')
                fstat.close()

            # se il test ha impiegato più di mezz'ora, da quello dopo stampiamo le statistiche
            if(resUsage["userTime"][0] >= statsThreshold and not printStats):
                print("starting to log statistics")
                printStats = True

        finally:
            newPoint['status'] = status
            # aggiungere tempo e memoria usati
            json.dump(newPoint, fout)
            fout.write(',\n')
            fout.close()





def main():
    #if not os.path.isfile("dataset.json"):
    #    fill()

    f = open("dataset.json", "r")
    #try:
    #    dataset = json.load(f)
    #except:
    #    f.close()
    #    fill()
    #    f = open("dataset.json", "r")
    #    dataset = json.load(f)

    dataset = json.load(f)
    #print(dataset[0])
    f.close()


    if(len(sys.argv) > 1):
        alg = sys.argv[1]       # è una stringa

        start, end = 0, len(dataset)
        if(len(sys.argv) > 2):
            start = int(sys.argv[2])
            if(len(sys.argv) > 3):
                end = int(sys.argv[3])

        testAlg(alg, dataset, start, end)

    else:
        for i in range(1, 5):
            testAlg(str(i), dataset)




if __name__=='__main__':
    main()