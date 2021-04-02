using Random, Distributions, CSV, DataFrames,StatsBase,Statistics,LinearAlgebra,DelimitedFiles

repy=parse(Int, ARGS[1]) 
Random.seed!(repy);  #different seed for different simulated phenotype

geno_path="/home/tianjing/paper_simu/pig.x.n928.p5024.txt"  #genotype of chromosome1 of pig
X=CSV.File(geno_path) |> DataFrame

SNP_ID= names(X)[2:end]
ID=X[:,1];

X=Matrix(X[:,2:end]);  #1st column is ID
n,p=size(X)


#genetic parameters
H2    = 0.5      # broad-sense heritability
var_y = 1        # phenotypic variance
var_e = (1-H2)*var_y
var_g = H2*var_y

mypi         = 0.95  # 1 - marker inclusion probability
nQTL         = floor(Int, p*(1-mypi))
QTL_position = sample(collect(1:p),nQTL);
QTL_ID       = SNP_ID[QTL_position]

open("QTL_position.txt", "w") do io
           writedlm(io, QTL_position, ',')
       end;
open("QTL_ID.txt", "w") do io
           writedlm(io, QTL_ID, ',')
       end;


# calculate Epistatic coefficient Interaction Ei
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

#calculate epistasis_coefficient*epistatic effect
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

Ma = X[:,QTL_position]

Md = copy(Ma)
Md[Md.!=1].=0.0;  # indicator for the heterozygous genotype

Mi = copy(Ma)

epistasis=create_epi(n, nQTL, ϵ, Mi)

#check h2
println("Vi:Vd:Va is",round(var(epistasis),digits=2),":",round(var(Md*d),digits=2),":",round(var(Ma*a),digits=2))

#scale Ma*a+Md*d+epistasis
g=Ma*a+Md*d+epistasis
println("Va+Vd+Vi before scaling:",round(var(g),digits=2))
g_new = g/std(g)*sqrt(var_g)
println("Va+Vd+Vi after scaling:",round(var(g_new),digits=2))

y= g_new +randn(n)*sqrt(var_e)

df = DataFrame(ID = ID, y = y);
CSV.write("pig.y.epistasis.rep$repy.txt", df)
