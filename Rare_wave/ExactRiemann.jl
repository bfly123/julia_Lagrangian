{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "Module ExactRiemann\n",
    "export Exact_Riemann\n",
    "function Presolve(varL::Var, varR::Var, conL::Const,conR::Const)\n",
    "        \n",
    "    ρL, uL, pL, sxxL = varL.ρ, varL.u, varL.p, varL.sxx\n",
    "    ρR, uR, pR, sxxR = varR.ρ, varR.u, varR.p, varR.sxx\n",
    "    σL  = -pL+sxxL\n",
    "    \n",
    "    UL = [ρL, uL, pL, sxxL]\n",
    "    UR = [ρR, uR, pR, sxxR]\n",
    "    cL = sound(UL[:],conL)\n",
    "   # eL = PToe(ρL, pL, conL)\n",
    "    \n",
    "    σR  = -pR+sxxR\n",
    "    cR = sound(UR[:],conR)\n",
    "  #  eR = PToe(ρR, pR, conR)\n",
    "    \n",
    "    sL = min(uL-cL, uR-cR)\n",
    "    sR = max(uL+cL, uR+cR)\n",
    "        \n",
    "    s_star = (σL - σR + ρL*uL*(sL-uL) - ρR* uR*(sR-uR))/(ρL*(sL-uL)-ρR*(sR-uR))\n",
    "        \n",
    "    ρLstar = ρL*(uL-sL)/(s_star - sL)\n",
    "    ρRstar = ρR*(uR-sR)/(s_star - sR)\n",
    "        \n",
    "    sxxLstar =  sxxL -4/3*conL.μ *(log(ρLstar) - log(ρL))\n",
    "    sxxRstar =  sxxR -4/3*conR.μ *(log(ρRstar) - log(ρR))\n",
    " #   @show sxxLstar, sxxRstar,ρLstar,ρRstar\n",
    "    caseL = CaseSelect(sxxL,sxxLstar, ρL, ρLstar,conL,1)\n",
    "    caseR = CaseSelect(sxxR,sxxRstar, ρR, ρRstar,conR,2)\n",
    "    σLstar = σL -ρL*(sL-uL)*(s_star-uL)\n",
    "    pLstar =sxxLstar - σLstar\n",
    "  #  @show pLstar,sxxLstar,ρLstar, s_star\n",
    "    \n",
    " #   return caseL, caseR\n",
    "    return ρLstar, ρRstar, caseL,caseR\n",
    "end\n",
    "\n",
    "function CaseSelect(sxx₀::Float64, sxx₁::Float64,ρ₀::Float64, ρ₁::Float64,con::Const,LoR::Int)\n",
    "    Y0 = con.Y0\n",
    "      if abs(sxx₁) ≥ 2/3*Y0  && abs(sxx₀) < 2/3*Y0   \n",
    "        if ρ₁ > ρ₀\n",
    "            if LoR == 1\n",
    "                case = Shock(2) # 左侧双激波\n",
    "            else\n",
    "                case = Shock(4) # 右侧双激波 \n",
    "            end\n",
    "        else\n",
    "            if LoR == 1\n",
    "                case = Rare(2) # case2 Left\n",
    "            else\n",
    "                case = Rare(4) # case2 Right\n",
    "            end\n",
    "        end\n",
    "    else\n",
    "        if ρ₁ > ρ₀\n",
    "            if LoR == 1\n",
    "                case = Shock(1) # case3\n",
    "                else\n",
    "                case = Shock(3)\n",
    "            end\n",
    "        else\n",
    "            if LoR == 1\n",
    "                case = Rare(1) #case1 Left\n",
    "            else\n",
    "                case = Rare(3)\n",
    "            end\n",
    "        end\n",
    "    end\n",
    "    return case\n",
    "end\n",
    "\n",
    "struct Shock\n",
    "    cs::Int\n",
    "end\n",
    "\n",
    "struct Rare\n",
    "    cs::Int\n",
    "    \n",
    "end\n",
    "##case 说明  cs=1 为左单波 cs=2 为左双波\n",
    "##          cs=3 为右单波 cs= 4 为右双波\n",
    "\n",
    "function Sxx(ρ, var0::Var, con::Const)\n",
    "    sxx0,ρ0 = var0.sxx,var0.ρ\n",
    "    Y0   = con.Y0\n",
    "    \n",
    "    sxx = sxx0 -4/3*con.μ *(log(ρ) - log(ρ0))\n",
    "    if abs(sxx) ≥ 2/3*Y0\n",
    "        sxx = 2/3*Y0*sign(sxx)\n",
    "    end\n",
    "    return sxx\n",
    "end\n",
    "\n",
    "function fp(ρ::Float64, var0::Var, con::Const, case::Rare)\n",
    "    ρ₀,p₀ = var0.ρ,var0.p\n",
    "    ρ0,Γ0,a0= con.ρ0,con.Γ0,con.a0   \n",
    "    λ₁ = ρ0*Γ0\n",
    "    f2(ρ) =a0^2*fηη(ρ, con) -λ₁*Sxx(ρ, var0,con)/ρ^2\n",
    "\n",
    "    p = p₀*exp(λ₁/ρ₀ - λ₁/ρ) + exp(-λ₁/ρ)*GaussIntegral(ρ->\n",
    "       f2(ρ)*exp(λ₁/ρ), ρ₀, ρ,5)\n",
    "    return p\n",
    "end\n",
    "#pLStar = fp(ρLStar,varL,con)\n",
    "\n",
    "function fu(ρ::Float64, var0::Var, con::Const,case::Rare)\n",
    "    ρ₀,u₀ = var0.ρ, var0.u\n",
    "    a0,ρ0,Γ0,μ= con.a0, con.ρ0, con.Γ0, con.μ\n",
    "    λ₁ = ρ0*Γ0\n",
    "    S2(ρ) = a0^2*fηη(ρ,con)+λ₁*(fp(ρ,var0,con,case)-Sxx(ρ,var0,con))/ρ^2\n",
    "    function f2(ρ)\n",
    "        if case.cs == 1 &&  case.cs ==2\n",
    "            return -√(S2(ρ)+4μ/(3ρ)) /ρ\n",
    "        else\n",
    "            return √(S2(ρ)+4μ/(3ρ)) /ρ\n",
    "        end\n",
    "    end\n",
    "    \n",
    "    u = u₀ + GaussIntegral(ρ->f2(ρ), ρ₀, ρ, 5)\n",
    "    return u\n",
    "end\n",
    "\n",
    "function GaussIntegral(f::Function,x₀::Float64,x₁::Float64,order::Int)\n",
    "    t₁= (x₁-x₀)/2\n",
    "    t₂= (x₁+x₀)/2\n",
    "    ω = zeros(Float64, 5)\n",
    "    p = zeros(Float64, 5)\n",
    "    \n",
    "    if order == 1\n",
    "        ω[1] = 2.0\n",
    "        p[1] = 0.0\n",
    "    elseif order == 3\n",
    "        ω[1] = 1.0; ω[2] = 1.0\n",
    "        p[1] = 1/√3.0; p[2] = -1/√3.0\n",
    "    elseif order == 5\n",
    "        ω[1] = 8.0/9; ω[2] = 5.0/9; ω[3] = 5.0/9\n",
    "        p[1] = 0.0; p[2] = -√(3.0/5); p[3] = √(3.0/5)\n",
    "    elseif order == 7\n",
    "        ω[1] = (18+√30)/36; ω[2] = (18+√30)/36\n",
    "        ω[3] = (18-√30)/36; ω[4] = (18-√30)/36\n",
    "        p[1] = √(3/7-2/7*√(6/5)); p[1] = -√(3/7-2/7*√(6/5))\n",
    "        p[3] = √(3/7+2/7*√(6/5)); p[4] = -√(3/7+2/7*√(6/5))\n",
    "    end\n",
    "    ∑ =sum( t₁*ω[i]*f(t₁*p[i]+t₂) for i in 1:floor(Int,order/2)+1)\n",
    "\n",
    "    return ∑\n",
    "end\n",
    "\n",
    "function Ṽar(var0::Var, case, con::Const)\n",
    "\n",
    "    Y0, μ= con.Y0, con.μ\n",
    "    sxx,ρ = var0.sxx,var0.ρ\n",
    "\n",
    "    cs = case.cs\n",
    "    if cs == 1 || cs == 3\n",
    "        var1 = var0\n",
    "    elseif cs == 2 || cs == 4\n",
    "        \n",
    "        if typeof(case) == Rare \n",
    "            sxx1 = 2/3*Y0\n",
    "            ρ1 = ρ*exp(-Y0/(2μ)+(3sxx)/(4μ))\n",
    "        elseif typeof(case) == Shock\n",
    "            sxx1 = -2/3*Y0\n",
    "            ρ1 = ρ*exp(Y0/(2μ)+(3sxx)/(4μ))\n",
    "        end\n",
    "        p1 = fp(ρ, var0, con, case)\n",
    "        u1 = fu(ρ, var0, con, case)\n",
    "    \n",
    "        var1= Var(ρ1, u1, p1, sxx1)\n",
    "    end\n",
    "    return var1\n",
    "end\n",
    "\n",
    "\n",
    "function fp(ρ::Float64, var0::Var, con::Const,case::Shock)\n",
    "   \n",
    "    ρ₀,p₀,sxx₀ = var0.ρ,var0.p,var0.sxx\n",
    "    ρ0,Γ0,a0 = con.ρ0,con.Γ0,con.a0\n",
    "    \n",
    "    σ₀ = -p₀+sxx₀\n",
    "    e₀ =  PToe(ρ₀, p₀, con)\n",
    "   # @show e₀\n",
    "    t =  (ρ-ρ₀)/(ρ₀*ρ) ; c₀ = 1/(Γ0*ρ0); c₁= a0^2/Γ0\n",
    "  \n",
    "    p = (2(c₁*fη(ρ,con) + e₀) -t*(σ₀+Sxx(ρ,var0,con)))/(-t+2c₀)\n",
    " \n",
    "    \n",
    " #   @show t, (2c₀-t),c₀\n",
    "    return p\n",
    "end\n",
    "\n",
    "function fu(ρ::Float64, var0::Var, con::Const,case::Shock)\n",
    "   \n",
    "    ρ₀,u₀,p₀,sxx₀ = var0.ρ,var0.u,var0.p,var0.sxx\n",
    "    σ₀ = -p₀+sxx₀\n",
    "    \n",
    "    sxx=Sxx(ρ,var0,con)\n",
    "    p  = fp(ρ, var0, con, case)\n",
    "    σ = -p+sxx\n",
    "    \n",
    "    t =  (ρ-ρ₀)/ρ₀*ρ\n",
    "    if case.cs == 1 || case.cs == 2\n",
    "        u = u₀ - √(t*(σ₀ - σ))\n",
    "    \n",
    "    else\n",
    "        u = u₀ + √(t*(σ₀ - σ))\n",
    "    end\n",
    "#    @show t, σ₀ - σ\n",
    "    return u\n",
    "end\n",
    "\n",
    "    fL(ρ) = fu(ρ, var0L, conL, caseL)\n",
    "    fR(ρ) = fu(ρ, var0R, conR, caseR)\n",
    "    gL(ρ) = -fp(ρ,var0L,conL,caseL)+Sxx( ρ, var0L, conL)\n",
    "    gR(ρ) = -fp(ρ,var0R,conR,caseR)+Sxx( ρ, var0R, conR)\n",
    "\n",
    "    fL¹(ρ) = Derivative(fL,ρ)\n",
    "    fR¹(ρ) = Derivative(fR,ρ)\n",
    "    gL¹(ρ) = Derivative(gL,ρ)\n",
    "    gR¹(ρ) = Derivative(gR,ρ)\n",
    "\n",
    "function Derivative(f::Function,x::Float64)\n",
    "    ϵ₀ = 1e-8\n",
    "    if abs(x) >= ϵ₀\n",
    "        ϵ= ϵ₀*x\n",
    "    else\n",
    "        ϵ =ϵ₀\n",
    "    end\n",
    "    f¹(x) = (f(x+ϵ)-f(x))/ϵ\n",
    "    return f¹(x)\n",
    "end\n",
    "\n",
    "function NewtonIter!(U::Array{Float64,1},J::Array{Float64,2},F::Array{Float64,1}, TOL::Float64 = 1.e-6) \n",
    "    I, = size(U)\n",
    "    U_new = zeros(Float64, I)\n",
    "    while true  \n",
    "        c = max(abs(F[1]),abs(F[2]))\n",
    "        if c <= TOL\n",
    "            break\n",
    "        end\n",
    "        U_new = U - inv(J)*F # same to inv(J)×F\n",
    "   #     c= CHA(U,U_new,F)\n",
    "     #   if c <= TOL\n",
    "      #      break\n",
    "      #  end\n",
    "    end\n",
    "    return U\n",
    "    \n",
    "end   \n",
    "\n",
    "function CHA(R0::Array{Float64,1}, R1::Array{Float64,1}, F::Array{Float64,1})\n",
    "    \n",
    "    ρₗ₀ = R0[1]; ρᵣ₀ = R0[2]\n",
    "    ρₗ₁ = R1[1]; ρᵣ₁ = R1[2]\n",
    "    \n",
    "    CHA = max(abs(2(ρₗ₁-ρₗ₀)/(ρₗ₁+ρₗ₀)),\n",
    "              abs(2(ρᵣ₁-ρᵣ₀)/(ρᵣ₁+ρᵣ₀)), F[1], F[2])\n",
    "    return CHA\n",
    "end   \n",
    "\n",
    "\n",
    "    fL(ρ, varL, conL, caseL) = fu(ρ, varL, conL, caseL)\n",
    "    fR(ρ, varR, conR, caseR) = fu(ρ, varR, conR, caseR)\n",
    "    gL(ρ, varL, conL, caseL) = -fp(ρ,varL,conL,caseL)+Sxx( ρ, varL, conL)\n",
    "    gR(ρ, varR, conR, caseR) = -fp(ρ,varR,conR,caseR)+Sxx( ρ, varR, conR)\n",
    "\n",
    "    fL¹(ρ, varL, conL, caseL) = Derivative(ρ->fL(ρ, varL, conL, caseL),ρ)\n",
    "    fR¹(ρ, varR, conR, caseR) = Derivative(ρ->fR(ρ, varR, conR, caseR),ρ)\n",
    "    gL¹(ρ, varL, conL, caseL) = Derivative(ρ->gL(ρ, varL, conL, caseL),ρ)\n",
    "    gR¹(ρ, varR, conR, caseR) = Derivative(ρ->gR(ρ, varR, conR, caseR),ρ)\n",
    " \n",
    "\n",
    "function Exact_Riemann(uL::Array{Float64,1}, uR::Array{Float64,1}, conL::Const, conR::Const )\n",
    "        \n",
    "    \n",
    "    var0L = Var(uL[1],uL[2],uL[3],uL[4])\n",
    "    var0R = Var(uR[1],uR[2],uR[3],uR[4])\n",
    " #  @show var0L\n",
    "    \n",
    "###预求解和分类    \n",
    "    ρLInit, ρRInit,caseL,caseR = Presolve(var0L, var0R, conL, conR)\n",
    "  # println(ρLInit,ρRInit)\n",
    "###求解  Q̂ \n",
    "    @show caseL,caseR\n",
    "    varL = Ṽar(var0L, caseL,conL)\n",
    "    varR = Ṽar(var0R, caseR,conR)\n",
    "    ρL, uL, pL, sxxL = varL.ρ, varL.u, varL.p, varL.sxx\n",
    "    ρR, uR, pR, sxxR = varR.ρ, varR.u, varR.p, varR.sxx\n",
    "\n",
    "   # @show varL, varR\n",
    "    ###迭代初始化    \n",
    "    ρLStar = ρLInit\n",
    "    ρRStar = ρRInit\n",
    "   \n",
    "    \n",
    "    TOL = 2.e-5\n",
    "    i = 1\n",
    "    U = zeros(Float64, 2)\n",
    "    \n",
    "    println(ρLStar,\"  \",ρRStar)\n",
    "    while i < 30\n",
    "        i += 1\n",
    "\n",
    "        J = [[-fR¹(ρRStar,varR, conR, caseR) fL¹(ρLStar,varL, conL, caseL)]; \n",
    "             [-gR¹(ρRStar,varR, conR, caseR) gL¹(ρLStar,varL, conL, caseL)]]\n",
    "    \n",
    "        U = [ρRStar, ρLStar]\n",
    "    \n",
    "        F = [fL(ρLStar,varL, conL, caseL)- fR(ρRStar,varR, conR, caseR), \n",
    "             gL(ρLStar,varL, conL, caseL)- gR(ρRStar,varR, conR, caseR)]\n",
    "      \n",
    "    \n",
    "      \n",
    "     c = max(abs(F[1]),abs(F[2])) \n",
    "    \n",
    "        @show U, F,J\n",
    "    \n",
    "      if c <= TOL\n",
    "            break\n",
    "        end\n",
    "    #   @show F\n",
    "        U = U - J\\F\n",
    "        ρRStar = U[1]\n",
    "        ρLStar = U[2]\n",
    "   end\n",
    "   \n",
    "    \n",
    "    ### 重新求解\n",
    "    pLStar = fp(ρLStar, varL, conL,caseL)\n",
    "    s_Star = fu(ρLStar, varL, conL,caseL)\n",
    "    sxxLStar = Sxx(ρLStar, varL, conL)\n",
    "\n",
    "    pRStar = fp(ρRStar, varR, conR,caseR)\n",
    "    sxxRStar = Sxx(ρRStar, varR, conR)\n",
    "    FL=zeros(Float64, 4)\n",
    "    FR=zeros(Float64, 4)\n",
    "    FL[1] = 0.0\n",
    "    FL[2] = pLStar-sxxLStar\n",
    "    FL[3] = (pLStar-sxxLStar)*s_Star\n",
    "    FL[4] = -4conL.μ/3*s_Star\n",
    "    uuh = s_Star\n",
    "\n",
    "    FR[1] = 0.0\n",
    "    FR[2] = pRStar-sxxRStar\n",
    "    FR[3] = (pRStar-sxxRStar)*s_Star\n",
    "    FR[4] = -4conR.μ/3*s_Star\n",
    "    uuh = s_Star\n",
    "    \n",
    "    return FL,FR, uuh\n",
    "end\n",
    "\n",
    "end #module"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Julia 1.1.0",
   "language": "julia",
   "name": "julia-1.1"
  },
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.1.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
