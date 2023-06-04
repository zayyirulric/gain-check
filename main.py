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
    # eq 1 test
    r_x = ((r_bp**(-1)) + (r_m**(-1)) + (r_mr**(-1)))**(-1)
    r_y = ((r_m**(-1)) + (r_bp**(-1)))**(-1)
    r_y = r_x

    num1 = (g_1+g_n+g_m)/(g_1+g_n+g_m+(1/r_1)+(1/r_y)+(1/r_n))
    num2 = (g_1+g_n)/(g_1+g_n+(1/r_1)+(1/r_n))
    den1 = ((1/r_bn)+(1/r_n)+(1/r_1))/(g_1+g_n+(1/r_n)+(1/r_1))
    den2 = ((1/r_n)+(1/r_1)+(1/r_y))/(g_1+g_n+g_m+(1/r_1)+(1/r_n)+(1/r_y))
    num = (num1-num2)
    den = den1-den2

    gain = num/den

    #################################### 
    # eq 2 test
    r_x = ((r_bp**(-1)) + (r_m**(-1)) + (r_mr**(-1)))**(-1)
    r_y = ((r_mr**(-1)) + (r_bp**(-1)))**(-1)

    n_1 = -((g_n/2)*(-g_n-(1/r_n)))/(g_n+(1/r_y)+(1/r_n))
    n_2 = -((g_1/2)*(g_1+(1/r_1)))/(g_1+g_m+(1/r_x)+(1/r_1))
    n_3 = -((g_m/2)*(g_1+(1/r_1)))/(g_1+g_m+(1/r_x)+(1/r_1))
    n_4 = g_1/2
    n_5 = -g_n/2

    d_1 = ((1/(2*r_n))*(-g_n-(1/r_n)))/(g_n+(1/r_y)+(1/r_n))
    d_2 = -((1/(2*r_1))*(g_1+(1/r_1)))/(g_1+g_m+(1/r_x)+(1/r_1))
    d_3 = 1/(2*r_bn)
    d_4 = 1/(2*r_n)
    d_5 = 1/(2*r_1)

    alt_num = n_1+n_2+n_3+n_4+n_5
    alt_den = d_1+d_2+d_3+d_4+d_5

    alt_gain = alt_num/alt_den

    print(f'{key} Gain (eqn. 1) = {gain}V/V, {db(gain)}dB')
    print(f'{key} Gain (eqn. 2) = {alt_gain}V/V, {db(alt_gain)}dB\n')