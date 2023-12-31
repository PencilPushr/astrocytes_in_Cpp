#include "probabilities.h"

double probabilities::phl[50] =
{
    0.333333333333, 0.0951756626745, 0.0271812917035, 0.00776661124715, 0.00222288659849, 0.000639921679991, 0.000187410416804, 5.81750667195e-05
, 2.17596143238e-05, 1.1436923158e-05, 7.88454209682e-06, 7.43738619183e-06, 7.37970057786e-06, 7.29603316347e-06, 7.27478942971e-06, 7.26006289992e-06
, 7.26084787208e-06, 7.26080132601e-06, 7.26061054996e-06, 7.26081361742e-06, 7.26079620991e-06, 7.26072365567e-06, 7.26058079345e-06, 7.26074725419e-06
, 7.26087576894e-06, 7.26073008288e-06, 7.26061194028e-06, 7.26074336727e-06, 7.26101528686e-06, 7.26081974085e-06, 7.26091667847e-06, 7.26058059059e-06
, 7.26084014577e-06, 7.2610063969e-06, 7.26069065682e-06, 7.26083092741e-06, 7.26076153595e-06, 7.26071756287e-06, 7.26092535023e-06, 7.26076421324e-06
, 7.26060026219e-06, 7.26075209967e-06, 7.26093367537e-06, 7.26073986493e-06, 7.26039032094e-06, 7.26091299989e-06, 7.26077756319e-06, 7.26071491915e-06
, 7.2607710224e-06, 7.26082337127e-06
};

// phl - probability of having low and high conducate - close - open state = closed gap junction
double probabilities::plh[50] =
{
    0.333333333333, 0.706523083189, 0.825908352326, 0.86117728508, 0.871356970354, 0.874272895353, 0.875107602476, 0.875346523333, 0.875414979813
, 0.875433066086, 0.875439533589, 0.875437420808, 0.875439605921, 0.875443525702, 0.875437685437, 0.875440602592, 0.875441218644, 0.875440938251
, 0.875438148928, 0.875441400815, 0.875441147657, 0.875440025483, 0.875437801975, 0.875440355641, 0.875442338151, 0.875440081279, 0.87543825197
, 0.875440285023, 0.87544449347, 0.875441466845, 0.875442967173, 0.875437765288, 0.875441782611, 0.875444355802, 0.875439468868, 0.875441639938
, 0.875440565916, 0.875439885313, 0.875443101387, 0.875440607355, 0.875438069767, 0.875440419864, 0.875443230241, 0.875440230498, 0.875434820356
, 0.875442910232, 0.875440813981, 0.875439844395, 0.875440712745, 0.875441522986
};

// phl - probability of having high and high conducate - open - open state = open gap junction
double probabilities::phh[50] =
{
    0.333333333333, 0.198301254137, 0.14691035597, 0.131056103673, 0.126420143048, 0.125087182967, 0.124704987107, 0.1245953016, 0.124563260573
, 0.124555496991, 0.124552581869, 0.124555141806, 0.124553014379, 0.124549178265, 0.124555039774, 0.124552137345, 0.124551520508, 0.124551800947
, 0.124554590462, 0.124551338371, 0.124551591547, 0.124552713793, 0.124554937445, 0.124552383612, 0.124550400973, 0.124552657991, 0.124554487418
, 0.124552454233, 0.124548245514, 0.124551272335, 0.12454977191, 0.124554974131, 0.124550956549, 0.124548383192, 0.124553270441, 0.124551099232
, 0.124552173322, 0.124552853969, 0.124549637687, 0.124552131881, 0.124554669633, 0.124552319384, 0.124549508825, 0.124552508762, 0.124557919254
, 0.124549828855, 0.124551925241, 0.12455289489, 0.124552026484, 0.124551216191
};