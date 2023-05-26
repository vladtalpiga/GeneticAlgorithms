import random
import math
import copy
from statistics import mean
import numpy as np
import matplotlib.pyplot as plt

# import struct
#
# def float_to_binary_string(value):
#     # Pack the float value as binary data
#     packed_value = struct.pack('!f', value)
#
#     # Convert the binary data to a binary string
#     binary_string = ''.join(format(byte, '08b') for byte in packed_value)
#
#     # Return the binary string truncated or padded to 22 characters
#     return binary_string[:22].ljust(22, '0')
#
# # Example usage
# float_value = 3.14159
# binary_string = float_to_binary_string(float_value)
# print(binary_string)


def dec_bin (nr, st, dr, p):
    l = math.ceil(math.log(((dr - st) * math.pow(10, p)), 2))
    d = (dr - st) / (math.pow(2, l))
    x = bin(math.floor((nr - a) / d))
    idx = x.find("b")
    print(idx)
    str = x[(idx+1):]
    x = str.rjust(22, "0")
    return x

def bin_dec (nr, st, dr, p):
    l = math.ceil(math.log(((dr - st) * math.pow(10, p)), 2))
    # d = (dr - st) / (math.pow(2, l)) +
    # nr = "0b" + nr
    # x = f"{(a + int(nr, 2) * d):.{l}f}"
    x = int(nr, 2)
    x = x * (dr-st)/(2 ** len(nr)) + st
    return float(x)

def calc_fitness(x, a, b, c):
    return float(a * math.pow(x, 2) + b * x + c)

def selectie(lfit):
    prob_selectie = []
    F = sum(lfit)
    s = 0
    intervale = [s/F]
    for i in range(len(lfit)):
        s += lfit[i]
        intervale.append(s/F)
        prob_selectie.append(lfit[i]/F)
    return intervale, prob_selectie

def binary_search(u, intervale, n):
    st = 0
    dr = n
    while st < dr:
        mij = (st + dr) // 2
        if intervale[mij] <= u < intervale[mij+1]:
            return mij + 1
        elif u < intervale[mij]:
            dr = mij
        else:
            st = mij + 1

def incrucisare(c1, c2, pt):
    return (c1[0][:pt] + c2[0][pt:]), (c2[0][:pt] + c1[0][pt:])

def mutatie(cromozomi, poz, prob):
    cr_bin = cromozomi[poz][0]
    k = False
    str = ""
    for i in range(len(cr_bin)):
        u = random.random()
        if u < prob:
            k = True
            if cr_bin[i] == '0':
                str += '1'
            else:
                str += '0'
        else:
            str += cr_bin[i]
    return str, k

def generat(n, l):
    pop = []
    for _ in range(n):
        bits = ''.join(random.choice(['0', '1']) for _ in range(l))
        pop.append(bits)
    return pop

### Date de intrare

f = open("input.txt", "r")
g = open("output.txt", "w")

n = int(f.readline())
st, dr = map(float, f.readline().split())
a, b, c = map(float, f.readline().split())
precizie = int(f.readline())
prob_recombinare = float(f.readline())
prob_mutatie = float(f.readline())
nr_etape = int(f.readline())

### Generam random populatia initiala de cromozomi

max_fitness = []
mean_fitness = []
cromozomi = []
fitness = []

g.write("Populatia initiala de cromozomi:\n")

l = math.ceil(math.log(((dr - st) * math.pow(10, precizie)), 2))
pop = generat(n, l)

i=0
for cr_bin in pop:
    cr_dec = bin_dec(cr_bin, st, dr, precizie)
    cromozomi.append((cr_bin, cr_dec))
    fit = calc_fitness(cr_dec, a, b, c)
    fitness.append(fit)
    g.write(f"{(5 - len(str(i + 1))) * ' ' + str(i + 1)}: {cr_bin}, x= {cr_dec}, f= {fit}\n")
    i +=1
# for b in pop:
#     cr_dec = random.uniform(st, dr)
#     cr_bin = dec_bin(cr_dec, st, dr, precizie)
#     print(cr_bin)
#     cromozomi.append((cr_bin, cr_dec))
#     fit = calc_fitness(cr_dec, a, b, c)
#     fitness.append(fit)



###

intervale_selectie, prob_selectie = selectie(fitness)

g.write("\nProbabilitati selectie cromozomi:\n")
i = 1
for ps in prob_selectie:
    g.write(f"{(5-len(str(i))) * ' ' + str(i)}: {ps}\n")
    i += 1

g.write("\nIntervale probabilitati selectie:\n")
i = 0
for ps in intervale_selectie:
    g.write(f"{(5-len(str(i))) * ' ' + str(i)}: {ps}\n")
    i += 1

g.write("\nEvidentierea procesului de selectie, folosind un numar aleator u uniform pe [0,1):\n")
cr_max = fitness.index(max(fitness)) + 1
cromozomi_selectati = [cr_max]
g.write(f"    1: Selectam cromozomul {cr_max}, care trece automat, deoarece are fitness-ul cel mai mare\n")
for i in range(1, n):
    u = random.random()
    cr_sel = binary_search(u, intervale_selectie, n)
    g.write(f"{(5-len(str(i+1))) * ' ' + str(i+1)}: u= {u} => Selectam cromozomul {cr_sel}\n")
    cromozomi_selectati.append(cr_sel)

cromozomi_dupa_selectie = []
fitness_dupa_selectie = []
g.write("\nCromozomii dupa selectie:\n")
for i in range(n):
    cr_sel = cromozomi_selectati[i]
    cr = cromozomi[cr_sel-1][0]
    cr_dec = cromozomi[cr_sel-1][1]
    fit = fitness[cr_sel-1]
    cromozomi_dupa_selectie.append((cr, cr_dec))
    fitness_dupa_selectie.append(fit)
    g.write(f"{(5 - len(str(i))) * ' ' + str(i)}: {cr}, x= {cr_dec}, f= {fit}\n")

cromozomi = copy.deepcopy(cromozomi_dupa_selectie)
fitness = copy.deepcopy(fitness_dupa_selectie)
participanti_recombinare = []

g.write(f"\nProbabilitatea de incrucisare: {prob_recombinare}")
for i in range(n):
    u = random.random()
    g.write(f"\n{(5 - len(str(i+1))) * ' ' + str(i+1)}: {cromozomi[i][0]}, x= {cromozomi[i][1]}, u= {u}")
    if u < prob_recombinare:
        participanti_recombinare.append(i+1)
        g.write(f" < {prob_recombinare} => Participa")

g.write("\n\n")
while len(participanti_recombinare) > 1:
    c1 = random.choice(participanti_recombinare)
    participanti_recombinare.remove(c1)
    c2 = random.choice(participanti_recombinare)
    participanti_recombinare.remove(c2)

    punct_taietura = random.randint(0, 22)

    g.write(f"Recombinare dintre cromozomul {c1} si cromozomul {c2}:\n   {cromozomi[c1-1][0]} x {cromozomi[c2-1][0]}, punct taietura: {punct_taietura}\n")
    aux1, aux2 = incrucisare(cromozomi[c1-1], cromozomi[c2-1], punct_taietura)
    cromozomi[c1-1] = (aux1, bin_dec(aux1, st, dr, precizie))
    cromozomi[c2-1] = (aux2, bin_dec(aux2, st, dr, precizie))
    g.write(f"=> {cromozomi[c1-1][0]} , {cromozomi[c2-1][0]}\n")

    fitness[c1-1] = calc_fitness(cromozomi[c1-1][1], a, b, c)
    fitness[c2-1] = calc_fitness(cromozomi[c2-1][1], a, b, c)

if len(participanti_recombinare) > 0:
    g.write(f"Cromozomul {participanti_recombinare[0]} nu a avut pereche pentru incrucisare.\n")

g.write("\nDupa recombinare:\n")
for i in range(n):
    g.write(f"{(5-len(str(i+1))) * ' ' + str(i+1)}: {cromozomi[i][0]}, x= {cromozomi[i][1]}, f= {fitness[i]}\n")

g.write(f"\nPosibilitate de mutatii pentru fiecare gena: {prob_mutatie}\n")

cromozomi_modificati = []
for i in range(n):
    cr_bin, k = mutatie(cromozomi, i, prob_mutatie)
    if k:
        cromozomi_modificati.append(i+1)
        cr_dec = bin_dec(cr_bin, st, dr, precizie)
        cromozomi[i] = (cr_bin, cr_dec)
        fitness[i] = calc_fitness(cr_dec, a, b, c)

g.write(f"Au fost modificati cromozomii: ")
for cr in cromozomi_modificati:
    g.write(f"{cr} ")

g.write("\n\nDupa mutatie:\n")
for i in range(n):
    g.write(f"{(5-len(str(i+1))) * ' ' + str(i+1)}: {cromozomi[i][0]}, x= {cromozomi[i][1]}, f= {fitness[i]}\n")

max_index = fitness.index(max(fitness))
max_fitness.append((cromozomi[max_index][1], fitness[max_index]))
mean_fitness.append(mean(fitness))

for j in range(nr_etape-1):
    intervale_selectie, prob_selectie = selectie(fitness)

    cr_max = fitness.index(max(fitness)) + 1
    cromozomi_selectati = [cr_max]
    for i in range(1, n):
        u = random.random()
        cr_sel = binary_search(u, intervale_selectie, n)
        cromozomi_selectati.append(cr_sel)

    cromozomi_dupa_selectie = []
    fitness_dupa_selectie = []
    for i in range(n):
        cr_sel = cromozomi_selectati[i]
        cr = cromozomi[cr_sel - 1][0]
        cr_dec = cromozomi[cr_sel - 1][1]
        fit = fitness[cr_sel - 1]
        cromozomi_dupa_selectie.append((cr, cr_dec))
        fitness_dupa_selectie.append(fit)

    cromozomi = copy.deepcopy(cromozomi_dupa_selectie)
    fitness = copy.deepcopy(fitness_dupa_selectie)
    participanti_recombinare = []

    for i in range(n):
        u = random.random()
        if u < prob_recombinare:
            participanti_recombinare.append(i + 1)

    while len(participanti_recombinare) > 1:
        c1 = random.choice(participanti_recombinare)
        participanti_recombinare.remove(c1)
        c2 = random.choice(participanti_recombinare)
        participanti_recombinare.remove(c2)

        punct_taietura = random.randint(0, 22)
        aux1, aux2 = incrucisare(cromozomi[c1 - 1], cromozomi[c2 - 1], punct_taietura)
        cromozomi[c1 - 1] = (aux1, bin_dec(aux1, st, dr, precizie))
        cromozomi[c2 - 1] = (aux2, bin_dec(aux2, st, dr, precizie))

        fitness[c1 - 1] = calc_fitness(cromozomi[c1 - 1][1], a, b, c)
        fitness[c2 - 1] = calc_fitness(cromozomi[c2 - 1][1], a, b, c)

    cromozomi_modificati = []
    for i in range(n):
        cr_bin, k = mutatie(cromozomi, i, prob_mutatie)
        if k:
            cromozomi_modificati.append(i + 1)
            cr_dec = bin_dec(cr_bin, st, dr, precizie)
            cromozomi[i] = (cr_bin, cr_dec)
            fitness[i] = calc_fitness(cr_dec, a, b, c)

    max_index = fitness.index(max(fitness))
    max_fitness.append((cromozomi[max_index][1], fitness[max_index]))
    mean_fitness.append(mean(fitness))

g.write("\n\nEvolutia maximului, respectiv a valorii medie:\n")

for i in range(nr_etape):
    g.write(f"{(5-len(str(i+1))) * ' ' + str(i+1)}: max= {max_fitness[i][1]}, mean= {mean_fitness[i]}\n")


### Plot

x = np.linspace(st, dr, int(dr - st) * 40)
y = np.vectorize(calc_fitness)(x, a, b, c)
plt.plot(x, y)

for x, y in max_fitness:
    plt.scatter(x, y)
    plt.pause(0.1)

plt.show()


f.close()
g.close()