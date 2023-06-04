import math

def metric_to_float(metric: str) -> float:
    metric_dict = {
        'p': 1e-12,
        'n': 1e-9,
        'u': 1e-6,
        'm': 1e-3,
        'k': 1e3,
        'M': 1e6,
        'G': 1e9
    }
    if metric[-1] in metric_dict:
        return float(metric[:-1]) * metric_dict[metric[-1]]
    else:
        return float(metric)
    
def db(gain: float) -> float:
    if (gain >= 0): pm = 1
    else: pm = -1
    db_gain = 20*math.log10(abs(gain))
    return pm*db_gain

param = {}

with open("params.txt","r") as f:
    key = ''
    for line in f.readlines():
        if (":" in line):
            key = line.split(":")[0]
            param[key] = {}
        else:        
            match line.split(" = ")[0]:
                case "g_1":
                    param[key]["g_1"] = metric_to_float(line.split(" = ")[1].strip())
                case "g_n":
                    param[key]["g_n"] = metric_to_float(line.split(" = ")[1].strip())
                case "g_m":
                    param[key]["g_m"] = metric_to_float(line.split(" = ")[1].strip())
                case "r_1":
                    param[key]["r_1"] = metric_to_float(line.split(" = ")[1].strip())
                case "r_n":
                    param[key]["r_n"] = metric_to_float(line.split(" = ")[1].strip())
                case "r_m":
                    param[key]["r_m"] = metric_to_float(line.split(" = ")[1].strip())
                case "r_bn":
                    param[key]["r_bn"] = metric_to_float(line.split(" = ")[1].strip())
                case "r_bp":
                    param[key]["r_bp"] = metric_to_float(line.split(" = ")[1].strip())
                case "r_mr":
                    param[key]["r_mr"] = metric_to_float(line.split(" = ")[1].strip())

for key in param.keys():
    g_1 = param[key]["g_1"]
    g_n = param[key]["g_n"]
    g_m = param[key]["g_m"]
    r_1 = param[key]["r_1"]
    r_n = param[key]["r_n"]
    r_m = param[key]["r_m"]
    r_bn = param[key]["r_bn"]
    r_bp = param[key]["r_bp"]
    r_mr = param[key]["r_mr"]

    ####################################
    # gain eq 1 test
    r_x = ((r_bp**(-1)) + (r_m**(-1)) + (r_mr**(-1)))**(-1)
    r_y = ((r_m**(-1)) + (r_bp**(-1)))**(-1)
    r_y = r_x

    n1_1 = (g_1+g_n+g_m)/(g_1+g_n+g_m+(1/r_1)+(1/r_y)+(1/r_n))
    n1_2 = (g_1+g_n)/(g_1+g_n+(1/r_1)+(1/r_n))
    d1_1 = ((1/r_bn)+(1/r_n)+(1/r_1))/(g_1+g_n+(1/r_n)+(1/r_1))
    d1_2 = ((1/r_n)+(1/r_1)+(1/r_y))/(g_1+g_n+g_m+(1/r_1)+(1/r_n)+(1/r_y))
    num1 = (n1_1-n1_2)
    den1 = d1_1-d1_2

    gain1 = num1/den1

    #################################### 
    # gain eq 2 test
    r_x = ((r_bp**(-1)) + (r_m**(-1)) + (r_mr**(-1)))**(-1)
    r_y = ((r_mr**(-1)) + (r_bp**(-1)))**(-1)

    n2_1 = -((g_n/2)*(-g_n-(1/r_n)))/(g_n+(1/r_y)+(1/r_n))
    n2_2 = -((g_1/2)*(g_1+(1/r_1)))/(g_1+g_m+(1/r_x)+(1/r_1))
    n2_3 = -((g_m/2)*(g_1+(1/r_1)))/(g_1+g_m+(1/r_x)+(1/r_1))
    n2_4 = g_1/2
    n2_5 = -g_n/2

    d2_1 = ((1/(2*r_n))*(-g_n-(1/r_n)))/(g_n+(1/r_y)+(1/r_n))
    d2_2 = -((1/(2*r_1))*(g_1+(1/r_1)))/(g_1+g_m+(1/r_x)+(1/r_1))
    d2_3 = 1/(2*r_bn)
    d2_4 = 1/(2*r_n)
    d2_5 = 1/(2*r_1)

    num2 = n2_1+n2_2+n2_3+n2_4+n2_5
    den2 = d2_1+d2_2+d2_3+d2_4+d2_5

    gain2 = num2/den2

    #################################### 
    # rout eq 1 test
    r_x = ((r_bp**(-1)) + (r_m**(-1)) + (r_mr**(-1)))**(-1)
    r_y = ((r_mr**(-1)) + (r_bp**(-1)))**(-1)

    d3_1 = ((1/r_1)/(g_1+g_m+(1/r_1)+(1/r_x)))*(g_1+(1/r_1))
    d3_2 = ((1/r_n)/(g_n+(1/r_n)+(1/r_y)))*(g_n+(1/r_n))
    
    den3 = (1/r_1) + (1/r_n) + (1/r_bn) - d3_1 - d3_2

    R_out1 = 1/den3

    gmro1 = -1*metric_to_float(key.split('S')[0])*R_out1

    #################################### 
    # rout eq 1 test
    r_x = ((r_bp**(-1)) + (r_m**(-1)) + (r_mr**(-1)))**(-1)
    r_y = ((r_mr**(-1)) + (r_bp**(-1)))**(-1)

    n4_1 = r_bn*r_1*r_n
    n4_2 = (g_1*r_1*r_x)+(g_m*r_1*r_x)+r_1+r_x
    n4_3 = (g_n*r_n*r_y)+r_n+r_y

    d4_1 = ((g_1*r_1*r_x)+(g_m*r_1*r_x)+r_1+r_x)*((g_n*r_n*r_y)+r_n+r_y)*((r_n*r_bn)+(r_1*r_bn)+(r_1*r_n))
    d4_2 = r_bn
    d4_3 = (r_n*r_x)*((g_1*r_1)+1)*((g_n*r_n*r_y)+r_n+r_y)
    d4_4 = (r_1*r_y)*((g_n*r_n)+1)*((g_1*r_1*r_x)+(g_m*r_1*r_x)+r_1+r_x)
    d4_5 = (-d4_3-d4_4)

    num4 = n4_1*n4_2*n4_3
    den4 = d4_1+(d4_2*d4_5)

    R_out2 = num4/den4

    gmro2 = -1*metric_to_float(key.split('S')[0])*R_out2

    print(f'{key} Gain (eqn. 1) = {gain1}V/V, {db(gain1)}dB')
    print(f'{key} Gain (eqn. 2) = {gain2}V/V, {db(gain2)}dB')
    print(f'{key} Rout (eqn. 1) = {R_out1}Ohms\n\tAv = -GmRout = {gmro1}V/V, {db(gmro1)}dB')
    print(f'{key} Rout (eqn. 2) = {R_out2}Ohms\n\tAv = -GmRout = {gmro2}V/V, {db(gmro2)}dB\n')