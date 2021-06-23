### Finding Normals Using RBFs and Tensors

## Including used packages
using LinearAlgebra
using SparseArrays
using NearestNeighbors
using DoubleFloats
using Plots
using LaTeXStrings
using BenchmarkTools
gr()


## Function definitions
# Circle
circX(t) = cos.(t);
circY(t) = sin.(t);

# True circle normals
circNormsX(t) = cos.(t);
circNormsY(t) = sin.(t);

# Pi curve
piX(t) = (17/31)*sin((235/57)+(-32)*t)+(19/17)*sin((192/55)+(-30)*t)+(47/32)*sin((69/25)+(-29)*t)+(35/26)*sin((75/34)+(-27)*t)+(6/31)*sin((23/10)+(-26)*t)+(35/43)*sin((10/33)+(-25)*t)+(126/43)*sin((421/158)+(-24)*t)+(143/57)*sin((35/22)+(-22)*t)+(106/27)*sin((84/29)+(-21)*t)+(88/25)*sin((23/27)+(-20)*t)+(74/27)*sin((53/22)+(-19)*t)+(44/53)*sin((117/25)+(-18)*t)+(126/25)*sin((88/49)+(-17)*t)+(79/11)*sin((43/26)+(-16)*t)+(43/12)*sin((41/17)+(-15)*t)+(47/27)*sin((244/81)+(-14)*t)+(8/5)*sin((79/19)+(-13)*t)+(373/46)*sin((109/38)+(-12)*t)+(1200/31)*sin((133/74)+(-11)*t)+(67/24)*sin((157/61)+(-10)*t)+(583/28)*sin((13/8)+(-8)*t)+(772/35)*sin((59/16)+(-7)*t)+(3705/46)*sin((117/50)+(-6)*t)+(862/13)*sin((19/8)+(-5)*t)+(6555/34)*sin((157/78)+(-3)*t)+(6949/13)*sin((83/27)+(-1)*t)+(-6805/54)*sin((1/145)+2*t)+(-5207/37)*sin((49/74)+4*t)+(-1811/58)*sin((55/43)+9*t)+(-63/20)*sin((2/23)+23*t)+(-266/177)*sin((13/18)+28*t)+(-2/21)*sin((7/16)+31*t);

piY(t) = (70/37)*sin((65/32)+(-32)*t)+(11/12)*sin((98/41)+(-31)*t)+(26/29)*sin((35/12)+(-30)*t)+(54/41)*sin((18/7)+(-29)*t)+(177/71)*sin((51/19)+(-27)*t)+(59/34)*sin((125/33)+(-26)*t)+(49/29)*sin((18/11)+(-25)*t)+(151/75)*sin((59/22)+(-24)*t)+(52/9)*sin((118/45)+(-22)*t)+(52/33)*sin((133/52)+(-21)*t)+(37/45)*sin((61/14)+(-20)*t)+(143/46)*sin((144/41)+(-19)*t)+(254/47)*sin((19/52)+(-18)*t)+(246/35)*sin((92/25)+(-17)*t)+(722/111)*sin((176/67)+(-16)*t)+(136/23)*sin((3/19)+(-15)*t)+(273/25)*sin((32/21)+(-13)*t)+(229/33)*sin((117/28)+(-12)*t)+(19/4)*sin((43/11)+(-11)*t)+(135/8)*sin((23/10)+(-10)*t)+(205/6)*sin((33/23)+(-8)*t)+(679/45)*sin((55/12)+(-7)*t)+(101/8)*sin((11/12)+(-6)*t)+(2760/59)*sin((40/11)+(-5)*t)+(1207/18)*sin((21/23)+(-4)*t)+(8566/27)*sin((39/28)+(-3)*t)+(12334/29)*sin((47/37)+(-2)*t)+(15410/39)*sin((185/41)+(-1)*t)+(-596/17)*sin((3/26)+9*t)+(-247/28)*sin((25/21)+14*t)+(-458/131)*sin((21/37)+23*t)+(-41/36)*sin((7/8)+28*t);

# True pi normals
piNormsX(t) = (abs((2240/37)*cos((65/32)+(-32)*t)+(341/12)*cos((98/41)+(-31)*t)+(780/29)*cos((35/12)+(-30)*t)+(1566/41)*cos((18/7)+(-29)*t)+(4779/71)*cos((51/19)+(-27)*t)+(767/17)*cos((125/33)+(-26)*t)+(1225/29)*cos((18/11)+(-25)*t)+(1208/25)*cos((59/22)+(-24)*t)+(1144/9)*cos((118/45)+(-22)*t)+(364/11)*cos((133/52)+(-21)*t)+(148/9)*cos((61/14)+(-20)*t)+(2717/46)*cos((144/41)+(-19)*t)+(4572/47)*cos((19/52)+(-18)*t)+(4182/35)*cos((92/25)+(-17)*t)+(11552/111)*cos((176/67)+(-16)*t)+(2040/23)*cos((3/19)+(-15)*t)+(3549/25)*cos((32/21)+(-13)*t)+(916/11)*cos((117/28)+(-12)*t)+(209/4)*cos((43/11)+(-11)*t)+(675/4)*cos((23/10)+(-10)*t)+(820/3)*cos((33/23)+(-8)*t)+(4753/45)*cos((55/12)+(-7)*t)+(303/4)*cos((11/12)+(-6)*t)+(13800/59)*cos((40/11)+(-5)*t)+(2414/9)*cos((21/23)+(-4)*t)+(8566/9)*cos((39/28)+(-3)*t)+(24668/29)*cos((47/37)+(-2)*t)+(15410/39)*cos((185/41)+(-1)*t)+(5364/17)*cos((3/26)+9*t)+(247/2)*cos((25/21)+14*t)+(10534/131)*cos((21/37)+23*t)+(287/9)*cos((7/8)+28*t))^2+abs((544/31)*cos((235/57)+(-32)*t)+(570/17)*cos((192/55)+(-30)*t)+(1363/32)*cos((69/25)+(-29)*t)+(945/26)*cos((75/34)+(-27)*t)+(156/31)*cos((23/10)+(-26)*t)+(875/43)*cos((10/33)+(-25)*t)+(3024/43)*cos((421/158)+(-24)*t)+(3146/57)*cos((35/22)+(-22)*t)+(742/9)*cos((84/29)+(-21)*t)+(352/5)*cos((23/27)+(-20)*t)+(1406/27)*cos((53/22)+(-19)*t)+(792/53)*cos((117/25)+(-18)*t)+(2142/25)*cos((88/49)+(-17)*t)+(1264/11)*cos((43/26)+(-16)*t)+(215/4)*cos((41/17)+(-15)*t)+(658/27)*cos((244/81)+(-14)*t)+(104/5)*cos((79/19)+(-13)*t)+(2238/23)*cos((109/38)+(-12)*t)+(13200/31)*cos((133/74)+(-11)*t)+(335/12)*cos((157/61)+(-10)*t)+(1166/7)*cos((13/8)+(-8)*t)+(772/5)*cos((59/16)+(-7)*t)+(11115/23)*cos((117/50)+(-6)*t)+(4310/13)*cos((19/8)+(-5)*t)+(19665/34)*cos((157/78)+(-3)*t)+(6949/13)*cos((83/27)+(-1)*t)+(6805/27)*cos((1/145)+2*t)+(20828/37)*cos((49/74)+4*t)+(16299/58)*cos((55/43)+9*t)+(1449/20)*cos((2/23)+23*t)+(7448/177)*cos((13/18)+28*t)+(62/21)*cos((7/16)+31*t))^2)^(-1/2)*((-2240/37)*cos((65/32)+(-32)*t)+(-341/12)*cos((98/41)+(-31)*t)+(-780/29)*cos((35/12)+(-30)*t)+(-1566/41)*cos((18/7)+(-29)*t)+(-4779/71)*cos((51/19)+(-27)*t)+(-767/17)*cos((125/33)+(-26)*t)+(-1225/29)*cos((18/11)+(-25)*t)+(-1208/25)*cos((59/22)+(-24)*t)+(-1144/9)*cos((118/45)+(-22)*t)+(-364/11)*cos((133/52)+(-21)*t)+(-148/9)*cos((61/14)+(-20)*t)+(-2717/46)*cos((144/41)+(-19)*t)+(-4572/47)*cos((19/52)+(-18)*t)+(-4182/35)*cos((92/25)+(-17)*t)+(-11552/111)*cos((176/67)+(-16)*t)+(-2040/23)*cos((3/19)+(-15)*t)+(-3549/25)*cos((32/21)+(-13)*t)+(-916/11)*cos((117/28)+(-12)*t)+(-209/4)*cos((43/11)+(-11)*t)+(-675/4)*cos((23/10)+(-10)*t)+(-820/3)*cos((33/23)+(-8)*t)+(-4753/45)*cos((55/12)+(-7)*t)+(-303/4)*cos((11/12)+(-6)*t)+(-13800/59)*cos((40/11)+(-5)*t)+(-2414/9)*cos((21/23)+(-4)*t)+(-8566/9)*cos((39/28)+(-3)*t)+(-24668/29)*cos((47/37)+(-2)*t)+(-15410/39)*cos((185/41)+(-1)*t)+(-5364/17)*cos((3/26)+9*t)+(-247/2)*cos((25/21)+14*t)+(-10534/131)*cos((21/37)+23*t)+(-287/9)*cos((7/8)+28*t));

piNormsY(t) = (abs((2240/37)*cos((65/32)+(-32)*t)+(341/12)*cos((98/41)+(-31)*t)+(780/29)*cos((35/12)+(-30)*t)+(1566/41)*cos((18/7)+(-29)*t)+(4779/71)*cos((51/19)+(-27)*t)+(767/17)*cos((125/33)+(-26)*t)+(1225/29)*cos((18/11)+(-25)*t)+(1208/25)*cos((59/22)+(-24)*t)+(1144/9)*cos((118/45)+(-22)*t)+(364/11)*cos((133/52)+(-21)*t)+(148/9)*cos((61/14)+(-20)*t)+(2717/46)*cos((144/41)+(-19)*t)+(4572/47)*cos((19/52)+(-18)*t)+(4182/35)*cos((92/25)+(-17)*t)+(11552/111)*cos((176/67)+(-16)*t)+(2040/23)*cos((3/19)+(-15)*t)+(3549/25)*cos((32/21)+(-13)*t)+(916/11)*cos((117/28)+(-12)*t)+(209/4)*cos((43/11)+(-11)*t)+(675/4)*cos((23/10)+(-10)*t)+(820/3)*cos((33/23)+(-8)*t)+(4753/45)*cos((55/12)+(-7)*t)+(303/4)*cos((11/12)+(-6)*t)+(13800/59)*cos((40/11)+(-5)*t)+(2414/9)*cos((21/23)+(-4)*t)+(8566/9)*cos((39/28)+(-3)*t)+(24668/29)*cos((47/37)+(-2)*t)+(15410/39)*cos((185/41)+(-1)*t)+(5364/17)*cos((3/26)+9*t)+(247/2)*cos((25/21)+14*t)+(10534/131)*cos((21/37)+23*t)+(287/9)*cos((7/8)+28*t))^2+abs((544/31)*cos((235/57)+(-32)*t)+(570/17)*cos((192/55)+(-30)*t)+(1363/32)*cos((69/25)+(-29)*t)+(945/26)*cos((75/34)+(-27)*t)+(156/31)*cos((23/10)+(-26)*t)+(875/43)*cos((10/33)+(-25)*t)+(3024/43)*cos((421/158)+(-24)*t)+(3146/57)*cos((35/22)+(-22)*t)+(742/9)*cos((84/29)+(-21)*t)+(352/5)*cos((23/27)+(-20)*t)+(1406/27)*cos((53/22)+(-19)*t)+(792/53)*cos((117/25)+(-18)*t)+(2142/25)*cos((88/49)+(-17)*t)+(1264/11)*cos((43/26)+(-16)*t)+(215/4)*cos((41/17)+(-15)*t)+(658/27)*cos((244/81)+(-14)*t)+(104/5)*cos((79/19)+(-13)*t)+(2238/23)*cos((109/38)+(-12)*t)+(13200/31)*cos((133/74)+(-11)*t)+(335/12)*cos((157/61)+(-10)*t)+(1166/7)*cos((13/8)+(-8)*t)+(772/5)*cos((59/16)+(-7)*t)+(11115/23)*cos((117/50)+(-6)*t)+(4310/13)*cos((19/8)+(-5)*t)+(19665/34)*cos((157/78)+(-3)*t)+(6949/13)*cos((83/27)+(-1)*t)+(6805/27)*cos((1/145)+2*t)+(20828/37)*cos((49/74)+4*t)+(16299/58)*cos((55/43)+9*t)+(1449/20)*cos((2/23)+23*t)+(7448/177)*cos((13/18)+28*t)+(62/21)*cos((7/16)+31*t))^2)^(-1/2)*((544/31)*cos((235/57)+(-32)*t)+(570/17)*cos((192/55)+(-30)*t)+(1363/32)*cos((69/25)+(-29)*t)+(945/26)*cos((75/34)+(-27)*t)+(156/31)*cos((23/10)+(-26)*t)+(875/43)*cos((10/33)+(-25)*t)+(3024/43)*cos((421/158)+(-24)*t)+(3146/57)*cos((35/22)+(-22)*t)+(742/9)*cos((84/29)+(-21)*t)+(352/5)*cos((23/27)+(-20)*t)+(1406/27)*cos((53/22)+(-19)*t)+(792/53)*cos((117/25)+(-18)*t)+(2142/25)*cos((88/49)+(-17)*t)+(1264/11)*cos((43/26)+(-16)*t)+(215/4)*cos((41/17)+(-15)*t)+(658/27)*cos((244/81)+(-14)*t)+(104/5)*cos((79/19)+(-13)*t)+(2238/23)*cos((109/38)+(-12)*t)+(13200/31)*cos((133/74)+(-11)*t)+(335/12)*cos((157/61)+(-10)*t)+(1166/7)*cos((13/8)+(-8)*t)+(772/5)*cos((59/16)+(-7)*t)+(11115/23)*cos((117/50)+(-6)*t)+(4310/13)*cos((19/8)+(-5)*t)+(19665/34)*cos((157/78)+(-3)*t)+(6949/13)*cos((83/27)+(-1)*t)+(6805/27)*cos((1/145)+2*t)+(20828/37)*cos((49/74)+4*t)+(16299/58)*cos((55/43)+9*t)+(1449/20)*cos((2/23)+23*t)+(7448/177)*cos((13/18)+28*t)+(62/21)*cos((7/16)+31*t));

# Circle function
circF(t) = sin(2*π*sin(t));

# True Laplace-Beltrami of Circle with f(t) = sin(2π sin(t))
trueCirc∇∇F(s) = ((-2*π*cos(2*π*sin(s))*sin(s))/sqrt(cos(s)^2 + sin(s)^2) - (4*π^2*cos(s)^2*sin(2*π*sin(s)))/sqrt(cos(s)^2 + sin(s)^2))/sqrt(cos(s)^2 + sin(s)^2);

# π function
piF(t) = sin(2*π*sin(t));

# True Laplace-Beltrami of π with f(t) = sin(2π sin(t))
truePi∇∇F(s) = (((-2240/37)*cos((65/32)+(-32)*s)+(-341/12)*cos((98/41)+(-31)*s)+(-780/29)*cos((35/12)+(-30)*s)+(-1566/41)*cos((18/7)+(-29)*s)+(-4779/71)*cos((51/19)+(-27)*s)+(-767/17)*cos((125/33)+(-26)*s)+(-1225/29)*cos((18/11)+(-25)*s)+(-1208/25)*cos((59/22)+(-24)*s)+(-1144/9)*cos((118/45)+(-22)*s)+(-364/11)*cos((133/52)+(-21)*s)+(-148/9)*cos((61/14)+(-20)*s)+(-2717/46)*cos((144/41)+(-19)*s)+(-4572/47)*cos((19/52)+(-18)*s)+(-4182/35)*cos((92/25)+(-17)*s)+(-11552/111)*cos((176/67)+(-16)*s)+(-2040/23)*cos((3/19)+(-15)*s)+(-3549/25)*cos((32/21)+(-13)*s)+(-916/11)*cos((117/28)+(-12)*s)+(-209/4)*cos((43/11)+(-11)*s)+(-675/4)*cos((23/10)+(-10)*s)+(-820/3)*cos((33/23)+(-8)*s)+(-4753/45)*cos((55/12)+(-7)*s)+(-303/4)*cos((11/12)+(-6)*s)+(-13800/59)*cos((40/11)+(-5)*s)+(-2414/9)*cos((21/23)+(-4)*s)+(-8566/9)*cos((39/28)+(-3)*s)+(-24668/29)*cos((47/37)+(-2)*s)+(-15410/39)*cos((185/41)+(-1)*s)+(-5364/17)*cos((3/26)+9*s)+(-247/2)*cos((25/21)+14*s)+(-10534/131)*cos((21/37)+23*s)+(-287/9)*cos((7/8)+28*s))^2+((-544/31)*cos((235/57)+(-32)*s)+(-570/17)*cos((192/55)+(-30)*s)+(-1363/32)*cos((69/25)+(-29)*s)+(-945/26)*cos((75/34)+(-27)*s)+(-156/31)*cos((23/10)+(-26)*s)+(-875/43)*cos((10/33)+(-25)*s)+(-3024/43)*cos((421/158)+(-24)*s)+(-3146/57)*cos((35/22)+(-22)*s)+(-742/9)*cos((84/29)+(-21)*s)+(-352/5)*cos((23/27)+(-20)*s)+(-1406/27)*cos((53/22)+(-19)*s)+(-792/53)*cos((117/25)+(-18)*s)+(-2142/25)*cos((88/49)+(-17)*s)+(-1264/11)*cos((43/26)+(-16)*s)+(-215/4)*cos((41/17)+(-15)*s)+(-658/27)*cos((244/81)+(-14)*s)+(-104/5)*cos((79/19)+(-13)*s)+(-2238/23)*cos((109/38)+(-12)*s)+(-13200/31)*cos((133/74)+(-11)*s)+(-335/12)*cos((157/61)+(-10)*s)+(-1166/7)*cos((13/8)+(-8)*s)+(-772/5)*cos((59/16)+(-7)*s)+(-11115/23)*cos((117/50)+(-6)*s)+(-4310/13)*cos((19/8)+(-5)*s)+(-19665/34)*cos((157/78)+(-3)*s)+(-6949/13)*cos((83/27)+(-1)*s)+(-6805/27)*cos((1/145)+2*s)+(-20828/37)*cos((49/74)+4*s)+(-16299/58)*cos((55/43)+9*s)+(-1449/20)*cos((2/23)+23*s)+(-7448/177)*cos((13/18)+28*s)+(-62/21)*cos((7/16)+31*s))^2)^(-1/2)*((-2)*pi*(((-2240/37)*cos((65/32)+(-32)*s)+(-341/12)*cos((98/41)+(-31)*s)+(-780/29)*cos((35/12)+(-30)*s)+(-1566/41)*cos((18/7)+(-29)*s)+(-4779/71)*cos((51/19)+(-27)*s)+(-767/17)*cos((125/33)+(-26)*s)+(-1225/29)*cos((18/11)+(-25)*s)+(-1208/25)*cos((59/22)+(-24)*s)+(-1144/9)*cos((118/45)+(-22)*s)+(-364/11)*cos((133/52)+(-21)*s)+(-148/9)*cos((61/14)+(-20)*s)+(-2717/46)*cos((144/41)+(-19)*s)+(-4572/47)*cos((19/52)+(-18)*s)+(-4182/35)*cos((92/25)+(-17)*s)+(-11552/111)*cos((176/67)+(-16)*s)+(-2040/23)*cos((3/19)+(-15)*s)+(-3549/25)*cos((32/21)+(-13)*s)+(-916/11)*cos((117/28)+(-12)*s)+(-209/4)*cos((43/11)+(-11)*s)+(-675/4)*cos((23/10)+(-10)*s)+(-820/3)*cos((33/23)+(-8)*s)+(-4753/45)*cos((55/12)+(-7)*s)+(-303/4)*cos((11/12)+(-6)*s)+(-13800/59)*cos((40/11)+(-5)*s)+(-2414/9)*cos((21/23)+(-4)*s)+(-8566/9)*cos((39/28)+(-3)*s)+(-24668/29)*cos((47/37)+(-2)*s)+(-15410/39)*cos((185/41)+(-1)*s)+(-5364/17)*cos((3/26)+9*s)+(-247/2)*cos((25/21)+14*s)+(-10534/131)*cos((21/37)+23*s)+(-287/9)*cos((7/8)+28*s))^2+((-544/31)*cos((235/57)+(-32)*s)+(-570/17)*cos((192/55)+(-30)*s)+(-1363/32)*cos((69/25)+(-29)*s)+(-945/26)*cos((75/34)+(-27)*s)+(-156/31)*cos((23/10)+(-26)*s)+(-875/43)*cos((10/33)+(-25)*s)+(-3024/43)*cos((421/158)+(-24)*s)+(-3146/57)*cos((35/22)+(-22)*s)+(-742/9)*cos((84/29)+(-21)*s)+(-352/5)*cos((23/27)+(-20)*s)+(-1406/27)*cos((53/22)+(-19)*s)+(-792/53)*cos((117/25)+(-18)*s)+(-2142/25)*cos((88/49)+(-17)*s)+(-1264/11)*cos((43/26)+(-16)*s)+(-215/4)*cos((41/17)+(-15)*s)+(-658/27)*cos((244/81)+(-14)*s)+(-104/5)*cos((79/19)+(-13)*s)+(-2238/23)*cos((109/38)+(-12)*s)+(-13200/31)*cos((133/74)+(-11)*s)+(-335/12)*cos((157/61)+(-10)*s)+(-1166/7)*cos((13/8)+(-8)*s)+(-772/5)*cos((59/16)+(-7)*s)+(-11115/23)*cos((117/50)+(-6)*s)+(-4310/13)*cos((19/8)+(-5)*s)+(-19665/34)*cos((157/78)+(-3)*s)+(-6949/13)*cos((83/27)+(-1)*s)+(-6805/27)*cos((1/145)+2*s)+(-20828/37)*cos((49/74)+4*s)+(-16299/58)*cos((55/43)+9*s)+(-1449/20)*cos((2/23)+23*s)+(-7448/177)*cos((13/18)+28*s)+(-62/21)*cos((7/16)+31*s))^2)^(-1/2)*cos(2*pi*sin(s))*sin(s)+(-1)*pi*cos(s)*(((-2240/37)*cos((65/32)+(-32)*s)+(-341/12)*cos((98/41)+(-31)*s)+(-780/29)*cos((35/12)+(-30)*s)+(-1566/41)*cos((18/7)+(-29)*s)+(-4779/71)*cos((51/19)+(-27)*s)+(-767/17)*cos((125/33)+(-26)*s)+(-1225/29)*cos((18/11)+(-25)*s)+(-1208/25)*cos((59/22)+(-24)*s)+(-1144/9)*cos((118/45)+(-22)*s)+(-364/11)*cos((133/52)+(-21)*s)+(-148/9)*cos((61/14)+(-20)*s)+(-2717/46)*cos((144/41)+(-19)*s)+(-4572/47)*cos((19/52)+(-18)*s)+(-4182/35)*cos((92/25)+(-17)*s)+(-11552/111)*cos((176/67)+(-16)*s)+(-2040/23)*cos((3/19)+(-15)*s)+(-3549/25)*cos((32/21)+(-13)*s)+(-916/11)*cos((117/28)+(-12)*s)+(-209/4)*cos((43/11)+(-11)*s)+(-675/4)*cos((23/10)+(-10)*s)+(-820/3)*cos((33/23)+(-8)*s)+(-4753/45)*cos((55/12)+(-7)*s)+(-303/4)*cos((11/12)+(-6)*s)+(-13800/59)*cos((40/11)+(-5)*s)+(-2414/9)*cos((21/23)+(-4)*s)+(-8566/9)*cos((39/28)+(-3)*s)+(-24668/29)*cos((47/37)+(-2)*s)+(-15410/39)*cos((185/41)+(-1)*s)+(-5364/17)*cos((3/26)+9*s)+(-247/2)*cos((25/21)+14*s)+(-10534/131)*cos((21/37)+23*s)+(-287/9)*cos((7/8)+28*s))^2+((-544/31)*cos((235/57)+(-32)*s)+(-570/17)*cos((192/55)+(-30)*s)+(-1363/32)*cos((69/25)+(-29)*s)+(-945/26)*cos((75/34)+(-27)*s)+(-156/31)*cos((23/10)+(-26)*s)+(-875/43)*cos((10/33)+(-25)*s)+(-3024/43)*cos((421/158)+(-24)*s)+(-3146/57)*cos((35/22)+(-22)*s)+(-742/9)*cos((84/29)+(-21)*s)+(-352/5)*cos((23/27)+(-20)*s)+(-1406/27)*cos((53/22)+(-19)*s)+(-792/53)*cos((117/25)+(-18)*s)+(-2142/25)*cos((88/49)+(-17)*s)+(-1264/11)*cos((43/26)+(-16)*s)+(-215/4)*cos((41/17)+(-15)*s)+(-658/27)*cos((244/81)+(-14)*s)+(-104/5)*cos((79/19)+(-13)*s)+(-2238/23)*cos((109/38)+(-12)*s)+(-13200/31)*cos((133/74)+(-11)*s)+(-335/12)*cos((157/61)+(-10)*s)+(-1166/7)*cos((13/8)+(-8)*s)+(-772/5)*cos((59/16)+(-7)*s)+(-11115/23)*cos((117/50)+(-6)*s)+(-4310/13)*cos((19/8)+(-5)*s)+(-19665/34)*cos((157/78)+(-3)*s)+(-6949/13)*cos((83/27)+(-1)*s)+(-6805/27)*cos((1/145)+2*s)+(-20828/37)*cos((49/74)+4*s)+(-16299/58)*cos((55/43)+9*s)+(-1449/20)*cos((2/23)+23*s)+(-7448/177)*cos((13/18)+28*s)+(-62/21)*cos((7/16)+31*s))^2)^(-3/2)*cos(2*pi*sin(s))*(2*((-2240/37)*cos((65/32)+(-32)*s)+(-341/12)*cos((98/41)+(-31)*s)+(-780/29)*cos((35/12)+(-30)*s)+(-1566/41)*cos((18/7)+(-29)*s)+(-4779/71)*cos((51/19)+(-27)*s)+(-767/17)*cos((125/33)+(-26)*s)+(-1225/29)*cos((18/11)+(-25)*s)+(-1208/25)*cos((59/22)+(-24)*s)+(-1144/9)*cos((118/45)+(-22)*s)+(-364/11)*cos((133/52)+(-21)*s)+(-148/9)*cos((61/14)+(-20)*s)+(-2717/46)*cos((144/41)+(-19)*s)+(-4572/47)*cos((19/52)+(-18)*s)+(-4182/35)*cos((92/25)+(-17)*s)+(-11552/111)*cos((176/67)+(-16)*s)+(-2040/23)*cos((3/19)+(-15)*s)+(-3549/25)*cos((32/21)+(-13)*s)+(-916/11)*cos((117/28)+(-12)*s)+(-209/4)*cos((43/11)+(-11)*s)+(-675/4)*cos((23/10)+(-10)*s)+(-820/3)*cos((33/23)+(-8)*s)+(-4753/45)*cos((55/12)+(-7)*s)+(-303/4)*cos((11/12)+(-6)*s)+(-13800/59)*cos((40/11)+(-5)*s)+(-2414/9)*cos((21/23)+(-4)*s)+(-8566/9)*cos((39/28)+(-3)*s)+(-24668/29)*cos((47/37)+(-2)*s)+(-15410/39)*cos((185/41)+(-1)*s)+(-5364/17)*cos((3/26)+9*s)+(-247/2)*cos((25/21)+14*s)+(-10534/131)*cos((21/37)+23*s)+(-287/9)*cos((7/8)+28*s))*((-71680/37)*sin((65/32)+(-32)*s)+(-10571/12)*sin((98/41)+(-31)*s)+(-23400/29)*sin((35/12)+(-30)*s)+(-45414/41)*sin((18/7)+(-29)*s)+(-129033/71)*sin((51/19)+(-27)*s)+(-19942/17)*sin((125/33)+(-26)*s)+(-30625/29)*sin((18/11)+(-25)*s)+(-28992/25)*sin((59/22)+(-24)*s)+(-25168/9)*sin((118/45)+(-22)*s)+(-7644/11)*sin((133/52)+(-21)*s)+(-2960/9)*sin((61/14)+(-20)*s)+(-51623/46)*sin((144/41)+(-19)*s)+(-82296/47)*sin((19/52)+(-18)*s)+(-71094/35)*sin((92/25)+(-17)*s)+(-184832/111)*sin((176/67)+(-16)*s)+(-30600/23)*sin((3/19)+(-15)*s)+(-46137/25)*sin((32/21)+(-13)*s)+(-10992/11)*sin((117/28)+(-12)*s)+(-2299/4)*sin((43/11)+(-11)*s)+(-3375/2)*sin((23/10)+(-10)*s)+(-6560/3)*sin((33/23)+(-8)*s)+(-33271/45)*sin((55/12)+(-7)*s)+(-909/2)*sin((11/12)+(-6)*s)+(-69000/59)*sin((40/11)+(-5)*s)+(-9656/9)*sin((21/23)+(-4)*s)+(-8566/3)*sin((39/28)+(-3)*s)+(-49336/29)*sin((47/37)+(-2)*s)+(-15410/39)*sin((185/41)+(-1)*s)+(48276/17)*sin((3/26)+9*s)+1729*sin((25/21)+14*s)+(242282/131)*sin((21/37)+23*s)+(8036/9)*sin((7/8)+28*s))+2*((-544/31)*cos((235/57)+(-32)*s)+(-570/17)*cos((192/55)+(-30)*s)+(-1363/32)*cos((69/25)+(-29)*s)+(-945/26)*cos((75/34)+(-27)*s)+(-156/31)*cos((23/10)+(-26)*s)+(-875/43)*cos((10/33)+(-25)*s)+(-3024/43)*cos((421/158)+(-24)*s)+(-3146/57)*cos((35/22)+(-22)*s)+(-742/9)*cos((84/29)+(-21)*s)+(-352/5)*cos((23/27)+(-20)*s)+(-1406/27)*cos((53/22)+(-19)*s)+(-792/53)*cos((117/25)+(-18)*s)+(-2142/25)*cos((88/49)+(-17)*s)+(-1264/11)*cos((43/26)+(-16)*s)+(-215/4)*cos((41/17)+(-15)*s)+(-658/27)*cos((244/81)+(-14)*s)+(-104/5)*cos((79/19)+(-13)*s)+(-2238/23)*cos((109/38)+(-12)*s)+(-13200/31)*cos((133/74)+(-11)*s)+(-335/12)*cos((157/61)+(-10)*s)+(-1166/7)*cos((13/8)+(-8)*s)+(-772/5)*cos((59/16)+(-7)*s)+(-11115/23)*cos((117/50)+(-6)*s)+(-4310/13)*cos((19/8)+(-5)*s)+(-19665/34)*cos((157/78)+(-3)*s)+(-6949/13)*cos((83/27)+(-1)*s)+(-6805/27)*cos((1/145)+2*s)+(-20828/37)*cos((49/74)+4*s)+(-16299/58)*cos((55/43)+9*s)+(-1449/20)*cos((2/23)+23*s)+(-7448/177)*cos((13/18)+28*s)+(-62/21)*cos((7/16)+31*s))*((-17408/31)*sin((235/57)+(-32)*s)+(-17100/17)*sin((192/55)+(-30)*s)+(-39527/32)*sin((69/25)+(-29)*s)+(-25515/26)*sin((75/34)+(-27)*s)+(-4056/31)*sin((23/10)+(-26)*s)+(-21875/43)*sin((10/33)+(-25)*s)+(-72576/43)*sin((421/158)+(-24)*s)+(-69212/57)*sin((35/22)+(-22)*s)+(-5194/3)*sin((84/29)+(-21)*s)+(-1408)*sin((23/27)+(-20)*s)+(-26714/27)*sin((53/22)+(-19)*s)+(-14256/53)*sin((117/25)+(-18)*s)+(-36414/25)*sin((88/49)+(-17)*s)+(-20224/11)*sin((43/26)+(-16)*s)+(-3225/4)*sin((41/17)+(-15)*s)+(-9212/27)*sin((244/81)+(-14)*s)+(-1352/5)*sin((79/19)+(-13)*s)+(-26856/23)*sin((109/38)+(-12)*s)+(-145200/31)*sin((133/74)+(-11)*s)+(-1675/6)*sin((157/61)+(-10)*s)+(-9328/7)*sin((13/8)+(-8)*s)+(-5404/5)*sin((59/16)+(-7)*s)+(-66690/23)*sin((117/50)+(-6)*s)+(-21550/13)*sin((19/8)+(-5)*s)+(-58995/34)*sin((157/78)+(-3)*s)+(-6949/13)*sin((83/27)+(-1)*s)+(13610/27)*sin((1/145)+2*s)+(83312/37)*sin((49/74)+4*s)+(146691/58)*sin((55/43)+9*s)+(33327/20)*sin((2/23)+23*s)+(208544/177)*sin((13/18)+28*s)+(1922/21)*sin((7/16)+31*s)))+(-4)*pi^2*cos(s)^2*(((-2240/37)*cos((65/32)+(-32)*s)+(-341/12)*cos((98/41)+(-31)*s)+(-780/29)*cos((35/12)+(-30)*s)+(-1566/41)*cos((18/7)+(-29)*s)+(-4779/71)*cos((51/19)+(-27)*s)+(-767/17)*cos((125/33)+(-26)*s)+(-1225/29)*cos((18/11)+(-25)*s)+(-1208/25)*cos((59/22)+(-24)*s)+(-1144/9)*cos((118/45)+(-22)*s)+(-364/11)*cos((133/52)+(-21)*s)+(-148/9)*cos((61/14)+(-20)*s)+(-2717/46)*cos((144/41)+(-19)*s)+(-4572/47)*cos((19/52)+(-18)*s)+(-4182/35)*cos((92/25)+(-17)*s)+(-11552/111)*cos((176/67)+(-16)*s)+(-2040/23)*cos((3/19)+(-15)*s)+(-3549/25)*cos((32/21)+(-13)*s)+(-916/11)*cos((117/28)+(-12)*s)+(-209/4)*cos((43/11)+(-11)*s)+(-675/4)*cos((23/10)+(-10)*s)+(-820/3)*cos((33/23)+(-8)*s)+(-4753/45)*cos((55/12)+(-7)*s)+(-303/4)*cos((11/12)+(-6)*s)+(-13800/59)*cos((40/11)+(-5)*s)+(-2414/9)*cos((21/23)+(-4)*s)+(-8566/9)*cos((39/28)+(-3)*s)+(-24668/29)*cos((47/37)+(-2)*s)+(-15410/39)*cos((185/41)+(-1)*s)+(-5364/17)*cos((3/26)+9*s)+(-247/2)*cos((25/21)+14*s)+(-10534/131)*cos((21/37)+23*s)+(-287/9)*cos((7/8)+28*s))^2+((-544/31)*cos((235/57)+(-32)*s)+(-570/17)*cos((192/55)+(-30)*s)+(-1363/32)*cos((69/25)+(-29)*s)+(-945/26)*cos((75/34)+(-27)*s)+(-156/31)*cos((23/10)+(-26)*s)+(-875/43)*cos((10/33)+(-25)*s)+(-3024/43)*cos((421/158)+(-24)*s)+(-3146/57)*cos((35/22)+(-22)*s)+(-742/9)*cos((84/29)+(-21)*s)+(-352/5)*cos((23/27)+(-20)*s)+(-1406/27)*cos((53/22)+(-19)*s)+(-792/53)*cos((117/25)+(-18)*s)+(-2142/25)*cos((88/49)+(-17)*s)+(-1264/11)*cos((43/26)+(-16)*s)+(-215/4)*cos((41/17)+(-15)*s)+(-658/27)*cos((244/81)+(-14)*s)+(-104/5)*cos((79/19)+(-13)*s)+(-2238/23)*cos((109/38)+(-12)*s)+(-13200/31)*cos((133/74)+(-11)*s)+(-335/12)*cos((157/61)+(-10)*s)+(-1166/7)*cos((13/8)+(-8)*s)+(-772/5)*cos((59/16)+(-7)*s)+(-11115/23)*cos((117/50)+(-6)*s)+(-4310/13)*cos((19/8)+(-5)*s)+(-19665/34)*cos((157/78)+(-3)*s)+(-6949/13)*cos((83/27)+(-1)*s)+(-6805/27)*cos((1/145)+2*s)+(-20828/37)*cos((49/74)+4*s)+(-16299/58)*cos((55/43)+9*s)+(-1449/20)*cos((2/23)+23*s)+(-7448/177)*cos((13/18)+28*s)+(-62/21)*cos((7/16)+31*s))^2)^(-1/2)*sin(2*pi*sin(s)));

## Generate node set given x and y coordinate sets
function dist(x, y)
    return [x'; y']
end

## Functions for generating approximate normals using a sorted node1
## distribution
# function for finding approximate normals
function approxNormals(nodes)
    # Number of nodes
    N = size(nodes,2);
    nmls = zeros(2,N);

    # First Node
    tmp = nodes[:,2] - nodes[:,N];
    tmp = rot90(tmp);
    nmls[:,1] = tmp;

    # Middle Nodes
    for i ∈ 2:N-1
        tmp = nodes[:,i+1] - nodes[:, i-1];
        tmp = rot90(tmp);
        nmls[:,i] = tmp;
    end

    # Last Node
    tmp = nodes[:,1] - nodes[:,N-1];
    tmp = rot90(tmp);
    nmls[:,N] = tmp;

    return nmls
end

# Function that rotates vectors -π/2
function rot90(vec)
    tmp = zeros(2);
    tmp[1] = vec[2];
    tmp[2] = -vec[1];

    return tmp
end

## Find indices for k nearest neighbors
# KNN for an unordered data set
function knnFull(nodes, n)
    kdtree = KDTree(nodes);
    idx, tmp = knn(kdtree, nodes, n, true);

    return idx
end

# KNN for ordered data set
function knnFullOrd(nodes, n)::Array{Array{Int,1},1}
    N = size(nodes,2);
    idx = [];
    for i ∈ 1:N
        push!(idx,zeros(n))
    end

    # Put nth node as nth first index
    for i ∈ 1:N
        idx[i][1] = i;
    end

    # Populate rest of indices
    for i ∈ 1:N, j ∈ 2:n
        if j%2 == 0
            k = i + floor(j/2);
            idx[i][j] = k>N ? k - N : k;
        else
            k = i - floor(j/2);
            idx[i][j] = k<=0 ? N+k : k;
        end
    end

    return idx
end


## Functions for setting up RBF envrionment
# Center main node
function center(nodes)
    nodes = nodes .- nodes[:,1];

    return nodes
end

# Finds angle between center vector and ĵ
sgn(x) = x >= 0 ? 1 : -1;
function findAngle(node)
    b = node[2]/norm(node);
    θ = sgn(node[1])acos(b);

    return θ
end

# Rotate nodes about origin
function rotate(nodes, θ)
    # Compute rotation matrix
    rotor = [cos(θ) -sin(θ);
             sin(θ) cos(θ)];

    # Rotate all vectors
    nodes = rotor*nodes;

    return nodes
end


## RBF related functions

# RBF definitons
# 1D
# Gaussian
# ϕ(x1,x2,m) = exp(-abs2(m*(x1-x2)));
# ϕ_x(x1,x2,m) = 2*(m^2)*(x2-x1)*exp(-abs2(m*(x1-x2)));

# PHS
ϕ(x1,x2,m) = abs(x1-x2)^m;
ϕ_x(x1,x2,m) = m*sign(x1-x2)*abs(x1-x2)^(m-1);
ϕ_xx(x1,x2,m) = m*(m-1)*abs(x1-x2)^(m-2);

# Interpolation weights for 1D with helping terms
function interpolate(nodes, m, o)
    n = size(nodes, 2);
    A0 = zeros(n,n);
    A1 = zeros(o+1,n);
    x = nodes[1,:];
    f = [nodes[2,:]; zeros(o+1)];

    # Define A without helping terms
    for i ∈ 1:n, j ∈ 1:n
        A0[i,j] = ϕ(x[i], x[j], m);
    end

    # Create helping terms
    A1[1,:] .= 1;
    for i ∈ 1:o, j ∈ 1:n
        A1[i+1,j] = x[j]^i;
    end

    # Create full A
    A = [A0 A1';
         A1 zeros(o+1,o+1)];

    # Solve for weights
    λ = A\f;

   return λ
end

# Interpolation function
function S(t, x, λ, m)
    n = size(x,1);
    nn = size(λ,1) - n;
    nt = size(t,1);

    s = zeros(nt);
    for i ∈ 1:n, j ∈ 1:nt
        s[j] += λ[i]*ϕ(t[j], x[i], m);
    end
    for i ∈ 0:nn-1, j ∈ 1:nt
        s[j] += λ[n+i+1]*t[j]^i;
    end

    return s
end

# RBF interpolant derivative at s=0
function S_x(x, λ, m)
    n = size(x,1);

    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_x(0, x[i], m);
    end

    # Add polynomial contribution if there is one
    if size(λ,1) > n+1
        s += λ[n+2];
    end

    return s
end

# RBF interpolant second derivative at s=0
function S_xx(x, λ, m)
    n = size(x,1);

    s = 0;
    for i ∈ 1:n
        s += λ[i]*ϕ_xx(0, x[i], m);
    end

    # Add polynomial contribution if there is one
    if size(λ,1) > n+2
        s += 2*λ[n+3];
    end

    return s
end

# Laplace-Beltrami Operator
function ∇∇(nodes, F, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);

    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFullOrd(nodes, n);

    # Allocate space for Laplace-Beltrami values
    vals = zeros(N);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights of local surface
        λs = interpolate(rot, m, o);

        # Compute RBF weights of local function
        λf = interpolate([rot[1,:]'; F[idx[i]]'], m, o);

        # Compute derivatives at S1 = 0
        S_s = S_x(rot[1,:], λs, m);
        S_ss = S_xx(rot[1,:], λs, m);
        S_f = S_x(rot[1,:], λf, m);
        S_ff = S_xx(rot[1,:], λf, m);

        # Compute length element
        s = 1 + S_s^2;

        # Compute ∇∇F
        vals[i] = s^(-1)*S_ff - s^(-2)*S_s*S_ss*S_f;
    end

    return vals
end

# Function for tensor normals
function findNormals(nodes, n, m, o)
    # Find number of nodes
    N = size(nodes, 2);

    # Compute approximate normals
    appNorms = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFullOrd(nodes, n);

    cent = zeros(2,n);
    rot = cent;
    normals = zeros(2,N);
    λ = zeros(n);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(appNorms[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights
        λ = interpolate(rot, m, o);

        # Compute derivative at S1 = 0
        f_s = S_x(rot[1,:], λ, m);

        # Compute length element
        s = 1 + f_s^2;
        L = sqrt(s);

        # Compute normal at S1 = 0
        tmp = [-f_s/L, 1/L];

        # Rotate local normal to proper angle
        normals[:,i] = rotate(tmp, -θ);
    end

    return normals
end

## Analysis functions
# Vector error
function vecError(trues, calcs)
    N = size(trues,2);
    normErrs = zeros(N);
    for j ∈ 1:N
        magDiff = norm(trues[:,j]-calcs[:,j]);
        mag = norm(trues[:,j]);
        normErrs[j] = magDiff/mag;
    end
    err = maximum(normErrs);
    return (err,normErrs)
end

# Scalar error
function scalError(trues, calcs)
    N = size(trues, 1);
    normErrs = zeros(N);
    for j ∈ 1:N
        magDiff = norm(trues[j]-calcs[j]);
        mag = norm(trues[j]);

        # A routine to handle mag ≈ 0
        if mag <= 10^(-13)
            normErrs[j] = 1magDiff;
        else
            normErrs[j] = magDiff/mag;
        end
    end
    err = maximum(normErrs);

    return (err,normErrs)
end

## Plot functions
# Plot of nodes and vectors
function vectorPlot(points, direction)
    N = size(points, 2);
    vecs = points + 75*direction;
    a = scatter(points[1,:],points[2,:],
                aspectratio = :equal,
                legend = false,
                show = false,
                markersize = 4,
                markerstrokewidth = 0,
                markerstrokealpha = 0,
                markeralpha = .9,
                title = L"\pi\mathrm{-Curve}\;\mathrm{Normals}");
    for i ∈ 1:N
        plot!([points[1,i], vecs[1,i]],
              [points[2,i], vecs[2,i]],
              color = :red,
              legend = false,
              show = false,
              lw = 1)
    end
    return a
end

# Nodes plot
function nodePlot(set1)
    a = scatter(set1[1,:],set1[2,:],
                color = :blue,
                show = false,
                legend = false,
                aspectratio = :equal,
                markersize = 1,
                markerstrokealpha = 0)
    return a
end

# Error plot
function errPlot(nodes, errs)
    N = size(errs,1);

    # Scale errors logarithmically
    errs = log10.(errs);

    a = scatter(nodes[1,:], nodes[2,:],
                marker_z = errs,
                c = :viridis,
                cbarlims = :auto,
                colorbar = :right,
                aspectratio = :equal,
                legend = false,
                markersize = 4,
                markerstrokealpha = 0.0,
                markerstrokewidth = 0.0,
                markeralpha = .75,
                title = L"\pi\mathrm{-Curve}\;\mathrm{Normal}\;\mathrm{Error}",
                colorbar_title = L"\log_{10}(Error)",
                dpi = 300)
    return a
end

# RBF interpolation plot
function interPlot(nodes, λ, m)
    x = nodes[1,:];
    y = nodes[2,:];
    t = range(minimum(x), maximum(x), length = 100);
    s = S(t, x, λ, m);

    a = scatter(x,y,
                color = :blue,
                aspectratio = :equal,
                show = false,
                legend = false,
                markersize = 1.5,
                markerstrokealpha = 0)
    plot!(t,s,lw=1)

    return a
end


# Animated process plot
function aniPlot(nodes, n = 10, m = 3, o = 3, delay = 0.001)
    # Find number of nodes
    N = size(nodes, 2);

    # Compute approximate normals
    nmls = approxNormals(nodes);

    # Find nearest neighbors
    idx = knnFullOrd(nodes, n);

    cent = zeros(2,n);
    rot = cent;
    λ = zeros(n);
    for i ∈ 1:N
        # Find angle
        θ = findAngle(nmls[:,i]);

        # Center nodes
        cent = center(nodes[:,idx[i]]);

        # Rotate nodes
        rot = rotate(cent, θ);

        # Compute RBF weights
        λ = interpolate(rot, m, o);

        # Plot data
        l = @layout [a b c; d e]
        a = plot(nodePlot(nodes),
                 title = "Initial Node Set",
                 xlims = (-1.1,1.1),
                 ylims = (-1.1,1.1))
        b = plot(nodePlot(nodes[:,idx[i]]),
                 title = "Nearest Neighbors",
                 xlims = (-1.1,1.1),
                 ylims = (-1.1,1.1))
        c = plot(nodePlot(cent),
                 title = "Centered",
                 xlims = (-.25,.25),
                 ylims = (-.25,.25))
        d = plot(nodePlot(rot),
                 title = "Rotated",
                 xlims = (-.25,.25),
                 ylims = (-.25,.25))
        e = plot(interPlot(rot, λ, m),
                 title = "Interpolated",
                 xlims = (-.1,.1),
                 ylims = (-.1,.1))

        f = plot(a,b,c,d,e,
                 layout = l)

        display(f)
        sleep(delay)
    end
end

## Main function for finding normals
function comp(N=100, n=10, m1=3, o=n-1)
    # Parameterizing our curve
    t = range(0,2*π-2*π/N, length = N);

    # Compute true normals
    truNorms = [piNormsX.(t)'; piNormsY.(t)'];
    for i ∈ 1:N
        truNorms[:,i] = truNorms[:,i]/norm(truNorms[:,i]);
    end

    # Generate nodes
    nodes = dist(piX.(t), piY.(t))
    appnorms = approxNormals(nodes)
    for i ∈ 1:N
        appnorms[:,i] = appnorms[:,i]/norm(appnorms[:,i]);
    end

    normals = findNormals(nodes, n, m1, o);

    # aniPlot(nodes, n, m1, o, 0.01)

    a = vectorPlot(nodes, normals);
    # b = vectorPlot(nodes, appnorms);
    # display(a)
    png(a, "Normals.png")

    a = vecError(truNorms, normals)[1]
    b = vecError(truNorms, appnorms)[1]

    c = errPlot(nodes, vecError(truNorms, normals)[2]);
    display(c)
    png(c, "Normal_Error.png")

    return [a,b]
end

function errs(m,o=-10)
    neighbors = [11];
    nodes = 1000:10000:100000;

    tmp = 0
    errapp = [];

    a = plot()
    a = plot(nodes,10^(-15.75)*nodes,
             label = "Rounding Error",
             linestyle = :dash,
             xaxis = :log,
             yaxis = :log,
             xlabel = "# Nodes",
             ylabel = "Inf-Norm Error",
             dpi = 300)
    a = plot!(nodes, ((9.7^3)./nodes).^11,
              label = latexstring("O(h^{11})"),
              linestyle = :dot,
              xaxis = :log,
              yaxis = :log)
    for n ∈ neighbors
        if o == -10
            oo = n-1;
        else
            oo = o;
        end

        err = [];
        for N ∈ nodes
            comps = comp(N, n, m, oo)
            push!(err,comps[1]);
            if tmp == 0
                push!(errapp,comps[2]);
            end
        end

        # if tmp == 0
        #     a = plot!(nodes,errapp,
        #               label = "Geometric",
        #               xaxis = :log,
        #               yaxis = :log)
        # end

        tmp = 1;

        a = plot!(nodes,err,
                  label = string(n," Neighbors"),
                  legend = :topright,
                  xaxis = :log,
                  yaxis = :log);
        display(a)
    end
end

# Laplace-Beltrami Stuff
function lapComp(N=100, n=10, m1=3, o=n-1)
    # Parameterizing our curve
    t = range(0,2*π-2*π/N, length = N);

    # Generate nodes and function values
    nodes = dist(piX.(t), piY.(t));
    F = piF.(t);

    # Compute true Laplace-Beltrami of F
    true∇∇F = truePi∇∇F.(t);

    laps = ∇∇(nodes, F, n, m1, o);

    a = scalError(true∇∇F, laps);

    # b = errPlot(nodes, a[2])
    # display(b)

    # c = plot3d(nodes[1,:],nodes[2,:],true∇∇F);
    # c = plot3d!(nodes[1,:],nodes[2,:],laps);
    # display(c)

    return a[1]
end

function lapErrs(m,o=-10)
    neighbors = [5,7,9,11,13];
    nodes = 1000:1000:20000;
    a = plot()
    a = plot(nodes,10^(-13.75)*nodes.^2,
             title = "RBF-Tensor Laplace-Beltrami(F)",
             label = "Rounding Error",
             linestyle = :dash,
             xaxis = :log,
             yaxis = :log,
             xlabel = "# Nodes",
             ylabel = "Inf-Norm Error",
             dpi = 300)
    # a = plot!(nodes, ((10^3)./nodes).^11.25,
    #           label = "Convergence Rate (L11.25)",
    #           xaxis = :log,
    #           yaxis = :log)
    for n ∈ neighbors
        if o == -10
            oo = n-1;
        else
            oo = o;
        end

        err = [];
        for N ∈ nodes
            comps = lapComp(N, n, m, oo)
            push!(err,comps);
        end

        a = plot!(nodes,err,
                  label = string(n, " Neighbors"),
                  legend = :topright,
                  xaxis = :log,
                  yaxis = :log);
        display(a)
    end
    # png(a, "fig_lap.png")
end

# Error plot
# errs(7)

# Single parameters
comp(1000,9,5,3)

# Laplace-Beltrami
# lapComp(1000, 11, 3, 5)

# lapErrs(7)
