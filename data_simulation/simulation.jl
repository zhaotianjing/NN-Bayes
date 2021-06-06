#this file is the Julia script to simulate data

using Random, Distributions, CSV, DataFrames,StatsBase,Statistics,LinearAlgebra,DelimitedFiles

### read data
repy=parse(Int, ARGS[1])   #e.g, y1
Random.seed!(repy);        #set different seed for y1,...,y20

geno_path="/home/tianjing/paper_simu/pig.x.n928.p5024.txt"  #path to genotype of chromosome1 of pig
X=CSV.File(geno_path) |> DataFrame                          #read genotype

SNP_ID= names(X)[2:end]   #get SNP ID (i.e., first row)
ID=X[:,1];                #get individual ID (i.e., first column)

X=Matrix(X[:,2:end]);     #get genotype matrix
n,p=size(X)               #n is #row, p is #column


### set genetic parameters
H2    = 0.5           # broad-sense heritability
var_y = 1             # phenotypic variance
var_e = (1-H2)*var_y  # environmental variance
var_g = H2*var_y      # genotipic variance

mypi         = 0.95                      # 5% markers are QTLs
nQTL         = floor(Int, p*(1-mypi))    #number of QTL
QTL_position = sample(collect(1:p),nQTL);#position of QTL
QTL_ID       = SNP_ID[QTL_position]      #corresponding SNP ID for QTL

### save data
open("QTL_position.txt", "w") do io
           writedlm(io, QTL_position, ',')
       end;
open("QTL_ID.txt", "w") do io
           writedlm(io, QTL_ID, ',')
       end;


#### calculate Epistatic coefficient (Ei)
# 
# van Bergen, Giel HH, et al. "Bayesian neural networks with variable selection 
# for prediction of genotypic values." Genetics Selection Evolution 52 (2020): 1-14.
#
function generate_epitasis_coef_Ec(qtl1,qtl2)
    n        = length(qtl1)
    epi_coef = zeros(n)  #length n
    for i in 1:n
        qtl_k = qtl1[i]
        qtl_l = qtl2[i]
        if qtl_k == 2            ##AA
            if qtl_l == 2      #BB
                epi_coef[i]=1
            elseif qtl_l == 1  #Bb
                epi_coef[i]=1
            else
                epi_coef[i]=0  #bb
            end
        elseif qtl_k == 1        ##Aa
            if qtl_l ==2      #BB
                epi_coef[i]=1
            elseif qtl_l==1    #Bb
                epi_coef[i]=1
            else
                epi_coef[i]=0 #bb
            end
        else                     ##aa
            if qtl_l ==2      #BB
                epi_coef[i]=0
            elseif qtl_l==1    #Bb
                epi_coef[i]=0
            else
                epi_coef[i]=0 #bb
            end
        end
    end
    return(epi_coef)
end

### calculate epistasis effect
function create_epi(n, nQTL, ϵ, Mi)
    epistasis=zeros(n)
    for i in 1:nQTL
        for j in i+1:nQTL
            ϵ_ij = ϵ[i,j] #epistatic effect for locus i,j
            qtli = Mi[:,i]
            qtlj = Mi[:,j]
            epistasis_coefficient_ij = generate_epitasis_coef_Ec(qtli,qtlj)
            epistasis +=  epistasis_coefficient_ij*ϵ_ij
        end
    end
    return epistasis
end


### simulation
# 
# van Bergen, Giel HH, et al. "Bayesian neural networks with variable selection 
# for prediction of genotypic values." Genetics Selection Evolution 52 (2020): 1-14.
#
a=randn(nQTL)                     #additive effect
δ=rand(Normal(1.2,0.3),nQTL)      #dominance factor
d=abs.(a).*δ                      #dominance effect

γ=triu(randn(nQTL,nQTL),1)        #epistatic factor

ϵ=zeros(nQTL,nQTL)                #epistatic effect
for i in 1:nQTL
    for j in i+1:nQTL
         ϵ[i,j]= γ[i,j]*sqrt(abs(a[i]*a[j]))
    end
end

Ma = X[:,QTL_position] #genotype matrix for QTL

Md = copy(Ma)
Md[Md.!=1].=0.0;       #indicator matrix for the heterozygous genotype

Mi = copy(Ma)
epistasis=create_epi(n, nQTL, ϵ, Mi)  #calculate epistasis


### scale Ma*a+Md*d+epistasis
g=Ma*a+Md*d+epistasis
println("Va+Vd+Vi before scaling:",round(var(g),digits=2))
g_new = g/std(g)*sqrt(var_g)
println("Va+Vd+Vi after scaling:",round(var(g_new),digits=2))

### simulate phenptype
y= g_new +randn(n)*sqrt(var_e)

### save data
df = DataFrame(ID = ID, y = y);
CSV.write("pig.y.epistasis.rep$repy.txt", df)
