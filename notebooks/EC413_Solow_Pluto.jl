### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 58753ecf-c322-4fb6-bb80-4abfc7cb5459
begin
	using PlutoUI
	#using PlotlyLight
	#using PyPlot
	using LaTeXStrings
	using Plots
	#pyplot()
	#Preset.Template.plotly_dark!()
	TableOfContents()
	#nothing
end

# ╔═╡ 1d082bd5-0086-4e29-aa5c-f4d41e1cf0fc
md"# _Solow Growth Model in Discrete Time (EC413-AT-2023)_

## How to use this notebook

This is an **interactive** _Julia_ notebook to let you play around with examples of Solow model that you have seen in the lecture and seminars.

#### Basic usage

- The good news is, you do not need to know _Julia_ (or any programming language) to make use of this notebook! For example, to see the effect of a change in parameters on the graphs or numerical results, simply use the sliders; changes will immediately take effect in the graphs and values (e.g., in the steady state level of capital per effective units of labour). If no slider is provided for the parameter you are interested in, you just need to change the parameter value at the point it is defined and then hit \"Shift+Enter\"

#### Advanced usage

- Feel free to edit this notebook as you wish! Change or delete the cells, values, functional forms etc., add new cells and experiment; break it and learn! **Don't worry, you won't be breaking the source code!** If you ended up breaking your copy too much, you can always close the browser and reload it again.
- In almost all cases, I make the cells (including either the markdown or th Julia code) invisible; If you're interested in reading, inspecting or even playing around with the code, hover your mouse to the left of a cell, there is an ``eye'' icon you can click on for the code to become visible; click it once more to make it invisible again.
- If you are interested in editing the code and you don't know Julia (but you know Matlab or Python), there is a shortcut: you can use the [excellent cheatsheet](https://cheatsheets.quantecon.org/) here. If you know R but not Julia, you can use this [R-Julia comparison cheatsheet](https://github.com/sswatson/cheatsheets/blob/master/jpr-cheatsheet.pdf).
- If you are going to inspect the code, notice that the Julia code in these notebooks are mainly for educational purposes, so it's not necessarily efficient or idiomatic. It's mostly written with the purpose that a wider range of the students can more easily understand the code.

Any questions, suggestions or issues? Don't hesitate to get in touch!

**Have fun!**
"

# ╔═╡ 0fee2ca0-7f57-4647-a187-36a9b3130cad
md"""
## The Baseline Model
Let's start by the Solow growth model. A Key building block of the Solow model is the aggregate production function, so first we define that. You can change it to experiment with other functional forms if you wish (again, press Shift+Enter for changes to take effect.)
"""

# ╔═╡ f8f03434-6a0a-4c91-a084-7f86a13eaf93
md"""
#### Aggregate Production Function
$\begin{align*}
Y_t&=  K_t^\alpha(A_t  N_t)^{1-\alpha}
\end{align*}$
"""

# ╔═╡ b2781e68-cc92-4186-a596-be753c52913e
F(A,K,N,α) = K^α*(A*N)^(1-α)

# ╔═╡ 3ff3caf6-a1e3-4298-ac99-2ff430443dc7
md"""
#### Laws of Motion
$\begin{align*}
K_{t+1}&=(1-\delta)K_t+s Y_t\\
N_{t+1}&=(1+n) N_t\\
A_{t+1}&= (1+g) A_t
\end{align*}$
"""

# ╔═╡ 9e596114-e8ef-40ca-aeca-12d02bc49f21
md"""
## BGP and the Steady State
### Transforming the Model and 
Letting $k_t=\frac{K_t}{A_t N_t}$ denote capital per effective unit of labour and combining the transformed law of motion of capital with the production function we get

$k_{t+1}=\frac{1}{(1+g)(1+n)}\left[(1-\delta)k_t +s  (k_t)^\alpha \right].$

Similarly, we can define output per effective unit of labour $y_t=\frac{Y_t}{A_t N_t}$, which is simply a function of $k_t$, $y_t = k_t^\alpha$.

**Q: why dividing by AN?**
- mathematics: to get a system which has a steady state (also dimension reduction)
- economics: Uzawa theorem
"""

# ╔═╡ f70f5669-04bf-495a-8adc-8bc7afab6e2e
md"""
### *The Uzawa Theorem (Optional)*
Assuming that
- aggregate production function $Y = F(K,N,A)$ is CRS (in $K, N$)
- Labour (population) $N$ grows with constant rate
- Resource constraint has the form $K_{t+1} = (1-\delta)K_t + Y_t - C_t$
then along a growth path where $Y$, $C$ and $K$ grow with *constant* positive rates,
1. we have $C$, $Y$ and $K$ growing with the **same** rate
2. technological progress can be represented in **purely labour augmenting form**

**Intuition** (and sketch of proof):
- *Asymmetry b/w K and N:* $K$ (but not $N$) accumulates using units of output $Y$; so $K$ *inherits* any growth in $Y$, grows with the same rate as $Y$, and $\frac{K_t}{Y_t}$ remains constant
- CRS implies $1=F(\frac{K_t}{Y_t}, \frac{N_t}{Y_t}, A_t)$ (balance production function),
- If $Y$ and $N$ grow with same rates $\frac{N}{Y}$ also remains constant, no technological progress
- but if $\frac{Y}{N}$ increases (as data shows), we need a compensating force for $N_t$ to keep up growing with $Y_t$ for the balance production function to be satisfied $\forall t$, i.e., purely labour-augmenting technological progress

Notes:
- This implies we can have non-Cobb-Douglas production functions giving us a BGP, provided that technological progress is solely labour-augmenting
- But in general and looking at the real world examples, no reason why all progress should be labor augmenting; in fact, many key innovations (steam engine; electricity; tractors; computers) are embodied in capital.
- That's why we might prefer to proceed with Cobb-Douglas that, thanks to having an elasticity of substitution between labour and capital equal to 1, is not only consistent with the stable factor shares but allows us to be agnostic about direction of technical change
"""

# ╔═╡ 66af3def-e62a-485a-ad93-901c3e920a96
md"""
### Parameters
Before going forward and to get a better sense of dynamics of the model and be able to do plotting and numerical experimentation, we need to parameterise the model. Let's start by the following:

$\begin{align*}
n &= 0\\
g &= 0.02\\
\delta &= 0.1\\
\alpha &= 0.3\\
s &= 0.15
\end{align*}$

- Note that in terms of data, s is between 10-25%, 
"""

# ╔═╡ 5669a9d4-c7cc-4426-a77d-33e23bc12175
begin
	n_slider = @bind n Slider(0:0.01:.1, default=0, show_value = true)
	g_slider = @bind g Slider(0:0.02:.2, default=.02, show_value = true)
	δ_slider = @bind δ Slider(0:0.05:1, default=.1, show_value = true)
	α_slider = @bind α Slider(0:0.05:1, default=.3, show_value = true)
	s_slider = @bind s Slider(0:0.05:1, default=.15, show_value = true)
	
	md"""
	**Golden-rule consumption as a function of α**
	
	n: $(n_slider)
	
	g: $(g_slider)
	
	δ: $(δ_slider)
	
	α: $(α_slider)
	
	s: $(s_slider)
	"""
end

# ╔═╡ cc25bdab-c8df-4640-ae4f-d58ba44fc263
md"""
### Finding the Steady State of the Transformed Model (analytically)
Steady state corresponds to $k_t$ not changing over time, i.e., $k_{t+1}=k_t = \bar{k}$; equation above implies:

$(g+n+gn+\delta)\bar{k} = s \bar{k}^\alpha$
$\bar{k}=\left(\frac{s}{g+n+gn+\delta}\right)^\frac{1}{1-\alpha}.$
"""

# ╔═╡ 29ee6c36-bcfb-40bc-bddf-7f8ba3d47c50
begin
	# finding the steady state value of ̄k in the transformed model
	k_ss = (s/(g+n+g*n+δ))^(1/(1-α))
	#println(k_ss)
	println("The steady state level of capital per effective units of labour is $(round(k_ss; digits=5))")
end

# ╔═╡ be361a68-cf13-4493-a810-f327c5f97f2c
md"""
### Finding the Steady State of the Transformed Model (visually)
There are at least two ways to inspect the steady state of this model using graphs.
1. As the intersection of the gross saving $sk^α$ and the depreciation line $(g+n+gn+\delta){k}$
2. As the intersection of $k_{t+1}$ (see it's law of motion above), and the 45ᵒ line.
"""

# ╔═╡ 3c9adefb-048c-4584-9ba8-b4ac22fdcd32
begin
	ks = 0:k_ss/32:k_ss*2
	save = k -> s*k^α
 	deprec= (n + g +n*g+δ)*ks
	plot(ks, save.(ks), label = "(gross) Saving", lw = 2, xlabel = L"$k$", ylabel = L"Saving/Depreciation",fmt=:png)
	plot!(ks, deprec, label = "Depreciation", ls=:dash, lw = 2)
	p_sav = plot!([k_ss], [save(k_ss)], seriestype = :scatter, label="Steady State")
	
	η = (1+g)*(1+n)
	ks_next = ((1-δ)*ks + save.(ks))/η
	plot(ks, ks_next, label = L"$k_{t+1}$", lw = 2, xlabel = L"$k_t$", ylabel = L"$k_{t+1}$",fmt=:png)
	p_next = plot!(ks, ks, label = "45ᵒ line", ls=:dash, lw = 2)
	p_next = plot!([k_ss], [k_ss], seriestype = :scatter, label="Steady State")

	plot(p_sav,p_next, layout=(1,2))
end

# ╔═╡ 53cc927f-02ec-4403-a3dc-6ea3791359b6
md"""
### Moments of Interest on the BGP
Now that we have found $k_t$ in the steady state, we can compute other moments of interest (which are a function of $k_t$ and parameters), evaluate them at the steady state of the transformed model, and go from that to the BGP of the original model.

- **Income per effective unit of labour:** $y_t = k_t^\alpha$; at the steady state/BGP: $\bar{y} = \bar{k}^\alpha$
- **Income per capita:** $\frac{Y_t}{N_t} = A_t y_t$; at the steady state/BGP, $\frac{Y_t}{N_t} = A_t \bar{y}$
- **Capital-output ratio:** $\frac{K_t}{Y_t} = \frac{k_t}{y_t} = k_t^{1-\alpha}$; at the steady state/BGP: $= \bar{k}^{1-\alpha}$
- **Rate of return to capital:** is equal to the marginal product of capital: $r_t = \frac{d Y_t}{d K_t} = \alpha k_t^{\alpha-1}$; at the steady state/BGP: $\bar{r} = \alpha \bar{k}^{\alpha-1}$
- **Capital Share:** $\frac{r_t K_t}{Y_t} = \alpha$ both at the transition path, and at the steady state/BGP
- **Labour share:** Similarly, $\frac{w_t N_t}{Y_t} = 1-\alpha$
"""

# ╔═╡ 66dbeb2a-32be-499a-a1d4-20c2647b8913
md"""
## Transition Path
Steady state is a point where if you start from, you will stay there forever. But what if you do not start from the steady state? I.e., what if a country starts from some initial level of capital $K_1$, productivity $A_1$, and population $N_1$ such that $k_1 = \frac{K_1}{A_1 N_1} \neq \bar{k}$?

Let's do a numerical exercise. Inspired by Assigment 2, let's consider an econmy, say England, which starts at period 1, with a low level of capital. Let's assume it starts with a level of capital per effective units of labour which is 1/1000 of it's steady state value.
"""

# ╔═╡ 2a5dd85f-8b53-4330-be69-3bee6b731caa
begin
	A₁ = 1;
	N₁ = 1;
	K₁ = A₁*N₁*(k_ss/1000)
	nothing
end

# ╔═╡ 76488fa6-2c94-4dbc-9c06-5bc4eaf6f3c6
md"""
What happens to the economy over time? Let's focus on the key variable of the (transformed) model, $k_t$, and follow its evolution for a horizon of 50. To this end, we first compute the output in the first period using the aggregate production function defined above:
"""

# ╔═╡ 3e2db772-2806-4470-b308-b6a8cd271e46
begin
	T = 50
	Y₁ = F(A₁, K₁, N₁, α)
	nothing
end

# ╔═╡ e417f1b3-fa3b-4736-80ef-d733fca46af8
md"""
Now we start from the initial values of A, K and N and use the laws of motion to compute future values of these. Moreover we can use the aggregate production function to compute the output in each period.

Here is a function to do this. You can ignore it unless you want to peek into the sausage factory!
"""

# ╔═╡ 97f40165-4cc3-4705-b19b-a73b5f10ff6d
function trajectories(A₁, K₁, N₁, α, g, n, δ, s, T)
#begin
	A = ones(T); A[1] = A₁ # defining a vector for productivities; initialising it with A₁
	N = ones(T); N[1] = N₁ # defining a vector for population; initialising it with N₁
	K = ones(T); K[1] = K₁ # defining a vector for capital; initialising it with K₁

	Y = ones(T); Y[1] = F(A₁, K₁, N₁, α) # defining a vector for output; initialising it with Y₁
	
	# now compute future values of variables of interest (A,N,K) according to their laws of motion
	for t ∈ 2:T
		A[t] = (1+g)*A[t-1]
		N[t] = (1+n)*N[t-1]
		K[t] = (1-δ)*K[t-1] + s*Y[t-1]
		Y[t] = F(A[t], K[t], N[t], α)
	end
	(A, K, N, Y)
end

# ╔═╡ c31be4ed-5d4d-4d9b-bd86-94c6020a3712
md"""
Now we can compute the full evolution of the economy (trajectories of Y, A, K and N, from period 1 to the last period) using the function defined above, and make use of it to further compute evolution of $k:$ 
"""

# ╔═╡ 02d9cee7-9ca7-48cf-8663-f7fbc23ebc00
md"""
What if we start from a higher level of capital? Although this is an empirically less relevant situation, let's say the initial level of capital is 10 times that of the steady state.
"""

# ╔═╡ 5d8e400e-9d1c-447c-8d4a-245d4f2a0fec
md"""
The graphs suggest that in both cases, $k_t$ converges back to its steady state value $\bar{k}$. This raises three questions:
1. is it a general property of the Solow model described above?
2. If yes, what is the reason?
3. And how fast such convergence happens. I.e., how each parameter of the model affects the speed at which convergence occurs?
"""

# ╔═╡ b4c56d87-3b1c-4a06-974b-87f733bc8018
md"""
## Convergence and Stability
### Convergence in the Solow Model
It is more straightforward to discuss convergence in terms of the transformed (rather than original) Solow model. Remember that the transformed Solow model is in effect a one-dimensional dynamical with a steady state. The criterion for having convergence in such a model is:

**If you start from the left of the steady state, you move to the right; if you start from the right, you move to the left.**
"""

# ╔═╡ 04a3af31-73a1-4de0-9f65-6ab2e4721447
md"""
To verify that this property holds for the Solow model, you can proceed either analytically or visually. For the former, one can compute $\Delta k_{t+1}$ and verify whether it is positive for $k_t < \bar{k}$ and negative for $k_t > \bar{k}$.

For the visual inspection, there are at least three ways to proceed, depicted in the graphs below:
"""

# ╔═╡ 9f6b6938-c148-48d3-ba78-39d6abd85403
begin
	#gr()
	ks_2 = 0:(k_ss/32):(k_ss*2)
	#save = k -> s*k^α
 	deprec_2= (n + g +n*g+δ)*ks_2
	plot(ks_2, save.(ks_2), label = "(gross) Saving", lw = 2, xlabel = L"$k$", ylabel = L"Saving/Depreciation", aspect_ratio = :none,fmt=:png)
	plot!(ks_2, deprec_2, label = "Depreciation", ls=:dash, lw = 2)
	p_sav_2 = plot!([k_ss], [save(k_ss)], seriestype = :scatter, label="Steady State")
	
	#η = (1+g)*(1+n)
	ks_next_2 = ((1-δ)*ks_2 + save.(ks_2))/η
	Δks = ks_next_2-ks_2
	plot(ks_2, ks_next_2, label = L"$k_{t+1}$", lw = 2, xlabel = L"$k_t$", ylabel = L"$k_{t+1}$",fmt=:png)
	p_next_2 = plot!(ks_2, ks_2, label = "45ᵒ line", ls=:dash, lw = 2)
	p_next_2 = plot!([k_ss], [k_ss], seriestype = :scatter, label="Steady State", aspect_ratio = :none)

	
	plot(ks_2, Δks, lw = 2, xlabel = L"$k_t$", ylabel = L"$\Delta k_{t+1}$", label = L"$\Delta k_{t+1}$",fmt=:png)
	plot!(ks_2,zeros(length(ks_2)), label = "zero", ls=:dash, lw = 2)
	p_dk = plot!([k_ss], [0], seriestype = :scatter, label="Steady State", aspect_ratio = :none)
	plot(p_sav_2,p_next_2, p_dk, layout=(2,2))
end

# ╔═╡ 30332b4c-1efb-4236-8284-076a36686164
md"""
1. Graph on the left: If we are to the left of the steady state (low capital), saving would be higher than depreciation, so capital (per effective units of labour) will increase and we will move to the right; vice versa if we start from a higher level of capital. $\implies$ convergence
2. Graph in the middle: if we start from the left of the steady state, next period capital per effective units of labour ($k_{t+1}$) is higher than the 45ᵒ line, which means $k_{t+1} > k_t$; so again we move to the right; and vice versa if we start from the right of the steady state. $\implies$ convergence
3. Graph on the right: If we start from the left of the steady state, $\Delta k_{t+1} = k_{t+1} - k_t > 0$, which amounts to the same thing, meaning we will move to the right and vice versa, $\implies$ convergence

Note: convergence and stability are very similar concepts which slightly differ in the usage: stability is usually with reference to the steady state, while convergence is with reference to an initial point which is usually distinct from the steady state.
"""

# ╔═╡ e4e7a385-d0f8-4648-a097-0c3b74c46e0c
md"""
### Speed of Convergence
Now that we showed the Solow model have the desired convergence property, we can proceed to the next question regarding the speed of convergence. There are several ways to quantify how fast the convergence occurs. Given our discrete time formulation, one way is to define the speed of convergence at period t, as the relative decline, in the relative gap between capital per effective unit of labour $k_t$ and the steady state $\bar{k}$ between periods $t$ and $(t+1)$:

$\frac{\frac{\bar{k} - k_t}{\bar{k}}- \frac{\bar{k} - k_{t+1}}{\bar{k}}}{\frac{\bar{k} - k_{t}}{\bar{k}}} = \frac{k_{t+1} - k_t}{{\bar{k}}-k_{t}}$

For example, in question 1 of Assignment 2 (where the steady state is modeled as another country, Scotland), there is a question about the speed of convergence for England in period one. In particular, assuming that England starts from a capital per effective unit of labour which is half that of the Scotland (the steady state), the speed of convergence between periods 1 and 2 will be equal to:

$\frac{\frac{1}{2}- \frac{k_{2,\text{Scotland}}-k_{2,\text{England}}}{k_{2,\text{Scotland}}}}{\frac{1}{2}}=1-2 \frac{k_{2,\text{Scotland}}-k_{2,\text{England}}}{k_{2,\text{Scotland}}}$
"""

# ╔═╡ 56f6341b-47df-4c86-b503-2999feb34ddc
begin
	n2_slider = @bind n2 Slider(0:0.01:.1, default=0, show_value = true)
	g2_slider = @bind g2 Slider(0:0.02:.2, default=.02, show_value = true)
	δ2_slider = @bind δ2 Slider(0:0.05:1, default=.1, show_value = true)
	α2_slider = @bind α2 Slider(0:0.05:1, default=.3, show_value = true)
	#s2_slider = @bind s2 Slider(0:0.05:1, default=.15, show_value = true)
	
	md"""
	**Speed of convergence between periods 1 and 2**
	
	n₂: $(n2_slider)
	
	g₂: $(g2_slider)
	
	δ₂: $(δ2_slider)
	
	α₂: $(α2_slider)
		
	"""
end

# ╔═╡ f48174cb-706a-4966-bf8f-20d07f288712
md"""
### Diminishing Marginal Product of Capital
Note that the key driver of convergence in the Solow model is the diminishing marginal product of capital. Higher levels of capital means lower return to capital and vice versa. This is important when you want to use the Solow model to think about transition, (conditional) convergence or episodes of rapid growth in the postwar economies (Japan, Germany, South Korea, etc.)
"""

# ╔═╡ a7beeea9-44a9-432b-b41c-32c67ba92afd
md"""
## The Kaldor Facts and the Solow Growth Model
### Facts of Growth
1. Near-constant growth in income per capita.
2. Capital-output ratio is nearly constant.
3. Return to capital is nearly constant.
4. Capital and labor share are nearly constant.

Note: These are facts of growth for developed/advanced economies.
Also note:
- Fact 1: note that income per capita grows with constant rate g in the BGP of the model; so growth in labour-augmenting productivity is the sole driver of growth in living standards in the baseline Solow model
- Fact 3: with competitive markets return to capital = MPK
"""

# ╔═╡ 6932442c-a3dd-4dfb-b5ed-1bc14ec39c7b
md"""
### The Solow Model and Facts of Growth
**Can the Solow model rationalise the facts of growth?**
Let's look at the moments of interest in the Solow model both along the transition path and the BGP?steady state:

| moment | Transition | BGP|
|:-------|:-----------|:---|
|growth rate of income per capita| $(1+g)\frac{k_{t+1}}{k_t}$| $1+g$|
|capital-output ratio | $k_t^{1-\alpha}$ | $\bar{k}^{1-\alpha}$|
|return to capital| $r_t = \alpha k_t^{\alpha-1}$| $\bar{r} = \alpha \bar{k}^{\alpha-1}$|
|labor share| $1-\alpha$|$1-\alpha$|
|capital share|$\alpha$|$\alpha$|

"""

# ╔═╡ 77287e52-6e9f-4228-8542-49fdffa9833f
md"""
- Apparently, facts of economic growth hold for the BGP of the model but not on the transition path (except for constancy of labour shares)
- This is consistent with the data: K/Y fluctuates quite a bit, but capital and labour share remain roughly constant
"""

# ╔═╡ 4ccdc633-e8eb-44bf-9d9f-95226b795ae3
md"""


### Why a Cobb-Douglas Production Function?
Why is the Aggregate Production Function Cobb-Douglas?

- This comes from Kaldor fact 4 (constant (non-degenerate) factor income shares), which calls for a unitary elasticity of substitution
- so 1% increase in interest rate relative to wages implies 1% decrease in capital-labour ratio, hence labour and capital share remain constant.



"""

# ╔═╡ 7fff17c8-b267-4d6f-b191-082ac3ac1a49
md"""
### Plotting Facts of Growth in the Solow Model
Again consider the economy we introduced above, which starts from a level of capital which is 1/1000 of the steady state. Let's plot the moments of interest regarding the growth facts.
"""

# ╔═╡ 5ccb172d-c2ee-4b94-b9d3-3b44b059dd1a
md"""
## The Golden Rule Consumption
*What is the highest consumption level that can be sustained in the Solow model?*

Note that the steady state income per effective units of labour is

$\bar{y}=\left(\frac{s}{g+n+gn+\delta}\right)^\frac{\alpha}{1-\alpha}.$

So on the BGP, income per capita is:

$\frac{Y_t}{N_t}=A_t \left(\frac{s}{g+n+gn+\delta}\right)^\frac{\alpha}{1-\alpha}.$

Saving rate is $s$, so consumption per capita is:

$\frac{C_t}{N_t}=A_t (1-s) \left(\frac{s}{g+n+gn+\delta}\right)^\frac{\alpha}{1-\alpha}.$

What saving rate $s^*$ maximises this expression? The expression on the RHS is strictly concave in s if $\alpha<\frac{1+s}{2}$ (empirically reasonable). The FOC implies that $s^* = \alpha$. So the Golden rule saving rate is equal to $\alpha$. Hence, the Golden rule level of capita per effective labour is:

${y^*}=\left(\frac{\alpha}{g+n+gn+\delta}\right)^\frac{\alpha}{1-\alpha}.$

Similarly the Golden rule consumption level of consumption per capita is:

$\frac{C_t^*}{N_t}=A_t (1-\alpha) \left(\frac{\alpha}{g+n+gn+\delta}\right)^\frac{\alpha}{1-\alpha}.$

The following graph shows how consumption per effective labour at BGP changes as a function of s (Assume $A_t = 1$):
"""

# ╔═╡ a486a61e-eb4d-4ad8-b7a6-fa3d95eeb682
begin
	α1_slider = @bind α1 Slider(0:0.05:1, default=.3, show_value = true)
	
	md"""
	**Golden-rule consumption as a function of α**
	
	α: $(α1_slider)
	
	"""
end

# ╔═╡ 0203b0b7-958a-42e3-a064-d601f165b5fc
begin
	grid_size = 100
	ss = range(0,1,length = grid_size)
	cs = zeros(grid_size)
	for i in eachindex(ss)
		si = ss[i]
		cs[i] = (1-si)*(si/(n+g+n*g+δ))^(α1/(1-α1))
	end
	
	sgold = α1
	c_golden = (1-sgold)*(sgold/(n+g+n*g+δ))^(α1/(1-α1))
	
	plot(ss, cs, xlabel = "s", ylabel="Consumption per capita", lw=2, label = "C/N",fmt=:png)
	#vline([c_golden])
	plot!([α1], seriestype = :vline, linestyle=:dash, label = "")
	plot!([α1],[c_golden], seriestype = :scatter, label="Golden-rule consumption")
end

# ╔═╡ 643c3a17-5536-4ed6-8faa-f5980ebd70e6
md"""
## Extensions
### Different Elasticities of Substitution (Between K, N)
We saw that, in the baseline Solow model, factor shares remain constant even on transition. This comes from a key propert of Cobb-Douglas production function: **unitary elasticity of substitution** between capital and labour.

What if we relax this assumption? The following production function generalises the Cobb-Douglas production function we have used thus far:

$Y = F(K,N, A) = \left( \alpha K^\rho + (1-\alpha) (AN)^\rho \right)^{\frac{1}{\rho}}$

With this production function, elasticity of substitution between K and N would be equal to

$\sigma = \frac{1}{1-\rho}$
That is, 1% increase in the the relative price of labour to capital ($\frac{w}{r}$) would translate to σ% increase in the capital labour ratio ($\frac{K}{N}$). Note that for $\rho=0$ (that is $\sigma=1$) we are back to the case of the Cobb-Douglas production function with capital elasticity $\alpha$ which we assumed at the start of the notebook. On the other hand, two extreme cases are:
- Linear production function has an elasticity of substitution equal to infinity;
- On the other hand, Leontieff production function (as in the Harrod-Domar model) has an elasticity of substitution equal to zero. Here, income share depends on which factor is abundant: if K abundant, then return to capital drops to zero, labour share becomes one; if (effective labour) is abundant, there would be mass unemployment, capital share is one.
"""

# ╔═╡ 1ec69369-e4a8-4e51-9925-bf502204350f
function trajectories(A₁, K₁, N₁, α, g, n, δ, s, σ, T, PF)
#begin
	A = ones(T); A[1] = A₁ # defining a vector for productivities; initialising it with A₁
	N = ones(T); N[1] = N₁ # defining a vector for population; initialising it with N₁
	K = ones(T); K[1] = K₁ # defining a vector for capital; initialising it with K₁

	Y = ones(T); Y[1] = PF(A₁, K₁, N₁, α, σ) # defining a vector for output; initialising it with Y₁
	
	# now compute future values of variables of interest (A,N,K) according to their laws of motion
	for t ∈ 2:T
		A[t] = (1+g)*A[t-1]
		N[t] = (1+n)*N[t-1]
		K[t] = (1-δ)*K[t-1] + s*Y[t-1]
		Y[t] = PF(A[t], K[t], N[t], α, σ)
	end
	(A, K, N, Y)
end

# ╔═╡ 031e0b88-732d-4ce2-9cd6-018faff82e84
begin
	A, K, N, Y = trajectories(A₁, K₁, N₁,α, g, n, δ, s, T)
	k = K./(A.*N)
	nothing
end

# ╔═╡ 0d6af4bc-011f-457a-b4d0-60eda96aa460
begin	
	plot(1:T, k, label = L"$k_t$ of England", lw = 2, xlabel = "time", ylabel = L"$k_t$",fmt=:png)
	plot!(1:T, k_ss*ones(T), label = L"Steady State $\bar{k}$", ls=:dash, lw = 2)
end

# ╔═╡ a279af3a-db5c-4f64-9036-a3679328896c
begin
	YN = Y./N
	YN_gr = YN[2:end]./YN[1:end-1] .- 1
	plot(1:T-1, YN_gr, label = "Over time", lw = 2, xlabel = L"$t$", ylabel = L"growth rate of $\frac{Y_{t}}{N_{t}}$", aspect_ratio = :none,fmt=:png)
	p_kaldor1 = plot!(1:T-1, g*ones(T-1), ls = :dash, lw=2, label = L"BGP")

	KY = K./Y
	plot(1:T-1, KY[1:T-1], label = "over time", lw = 2, xlabel = L"$t$", ylabel = L"$\frac{K_{t}}{Y_{t}}$", aspect_ratio = :none,fmt=:png)
	p_kaldor2 = plot!(1:T-1, (k_ss^(1-α))*ones(T-1), ls = :dash, lw=2, label = L"BGP")

	r = α*k.^(α-1)
	plot(1:T-1, r[1:T-1], label = L"$r_t$ (over time)", lw = 2, xlabel = L"$t$", ylabel = L"${r_t}$", aspect_ratio = :none,fmt=:png)
	p_kaldor3 = plot!(1:T-1, (α*k_ss^(α-1))*ones(T-1), ls = :dash, lw=2, label = L"$\bar{r}$ (BGP)")

	KS = r.*K./Y
	plot(1:T-1, KS[1:T-1], label = "Capital Share over time", lw = 2, xlabel = L"$t$", ylabel = L"$\frac{r_t K_t}{Y_t}$", aspect_ratio = :none, ylimits=(0,1),fmt=:png)
	p_kaldor4 = plot!(1:T-1, α*ones(T-1), ls = :dash, lw=2, label = "Capital share (BGP)")

	plot(p_kaldor1,p_kaldor2, p_kaldor3, p_kaldor4, layout=(2,2))
end

# ╔═╡ 17307ffa-4f1c-46e1-bb77-42c978464b91
begin
	K₁ₕ = A₁*N₁*(k_ss*10)
	Aₕ, Kₕ, Nₕ, Yₕ = trajectories(A₁, K₁ₕ, N₁,α, g, n, δ, s, T)
	kₕ = Kₕ./(Aₕ.*Nₕ)
	begin	
	plot(1:T, kₕ, label = L"$k_t$ of England", lw = 2, xlabel = "time", ylabel = L"$k_t$",fmt=:png)
	plot!(1:T, k_ss*ones(T), label = L"Steady State $\bar{k}$", ls=:dash, lw = 2)
	end
end

# ╔═╡ 900c31c6-12eb-4f4d-95c8-9aa4ec227beb
begin
	K₁ₑ = A₁*N₁*k_ss/2
	Aₑ, Kₑ, Nₑ, Yₑ = trajectories(A₁, K₁ₑ, N₁,α2, g2, n2, δ2, s, T)
	kₑ = Kₑ./(A.*N)
	
	speed_12 = (2^(1-α2)-1)*(1+(δ2-1)/((1+n2)*(1+g2)))
	println("The speed of convergence is $(round(speed_12; digits=5))")
end

# ╔═╡ 3cdd278a-1cb6-4ab7-92fb-7798a7eab02f
md"""
#### Elasticity of Substitution and Facts of Growth
Let's start again from a very low level of capital, 1/1000 of the steady state, and use the same parameters as before (except for the elasticity of substitution which you can change).

Now you can play around with the σ parameter and observe how moments of interest change both along the transition path and the BGP. In particular focus on how things change when you increase σ above one (higher substitutability) or towards zero (higher complementarity):
- How does the capital share change?
- What about the return to capital?

Another interesting set of observations relate to the transitional dynamics. In particular, decrease σ towards zero and observe how the speed of transition changes. Increase the T₂ slider to see whether and when (how fast) the transition happens. How do you inetrpret this?

(if you want to change other parameters in addition to σ, you should go back to the Parameters section above and use the sliders there)
"""

# ╔═╡ 9620ee95-69ca-44e3-892b-394b68a00a19
begin
	σ_slider = @bind σ Slider(0.:.2:10, default=1, show_value = true)
	T_slider = @bind T₂ Slider(10:10:500, default=50, show_value = true)
	md"""
	**Elasticity of substitution between K and N**
	
	σ: $(σ_slider)
	
	T₂: $(T_slider)  (number of transition periods)
	"""
end

# ╔═╡ 00650bbb-f830-45e9-ab38-00c034ac8c2b
function F_ces(A,K,N,α, σ)
	if σ==0.
		return min(K, A*N)
	elseif σ==1.
		return F(A,K,N,α)
	else
		ρ = 1-(1/σ)
		(α*K^ρ+(1-α)*((A*N)^ρ))^(1/ρ)
	end
end

# ╔═╡ ea573837-22e3-4f85-8ab6-08cd048dc1b1
begin
	A_ces, K_ces, N_ces, Y_ces = trajectories(A₁, K₁, N₁,α, g, n, δ, s, σ, T₂, F_ces)
	k_ces = K_ces./(A_ces.*N_ces)
	ρ = 1-(1/σ)
	k_ss_ces = σ==1 ? k_ss : s*((1-α)/((g+n+g*n+δ)^ρ - α*s^ρ))^(1/ρ)
	y_ss_ces = σ==1 ? k_ss^α : (α*k_ss_ces^(ρ)+(1-α))^(1/ρ)
	
	YN_ces = Y_ces./N_ces

	YN_gr_ces = YN_ces[2:end]./YN_ces[1:end-1] .- 1
	plot(1:T₂-1, YN_gr_ces, label = "Over time", lw = 2, xlabel = L"$t$", ylabel = L"growth rate of $\frac{Y_{t}}{N_{t}}$", aspect_ratio = :none)
	p_kaldor1_ces = plot!(1:T₂-1, g*ones(T₂-1), ls = :dash, lw=2, label = L"BGP")

	KY_ces = K_ces./Y_ces
	plot(1:T₂-1, KY_ces[1:T₂-1], label = "over time", lw = 2, xlabel = L"$t$", ylabel = L"$\frac{K_{t}}{Y_{t}}$", aspect_ratio = :none)
	p_kaldor2_ces = plot!(1:T₂-1, (k_ss_ces/y_ss_ces)*ones(T₂-1), ls = :dash, lw=2, label = L"BGP")

	if σ==0
		r_ces = zeros(T₂)
		for t=1:T₂
			r_ces[t] = K_ces[t] <= A_ces[t]*N_ces[t] ? 1 : 0
		end
		r_ss_ces = 0
	else
		r_ces = σ==1 ? α*k_ces.^(α-1) : α*(α.+(1-α)*k_ces.^(1/σ - 1)).^(1/(1-σ))
		r_ss_ces = σ==1 ? α*k_ss_ces^(α-1) : α*(α+(1-α)*k_ss_ces^(1/σ - 1))^(1/(1-σ))
	end	
	
	plot(1:T₂-1, r_ces[1:T₂-1], label = L"$r_t$ (over time)", lw = 2, xlabel = L"$t$", ylabel = L"${r_t}$", aspect_ratio = :none)
	p_kaldor3_ces = plot!(1:T₂-1, (r_ss_ces)*ones(T₂-1), ls = :dash, lw=2, label = L"$\bar{r}$ (BGP)")

	KS_ces = r_ces.*K_ces./Y_ces
	KS_ss_ces = r_ss_ces*k_ss_ces/y_ss_ces
	plot(1:T₂-1, KS_ces[1:T₂-1], label = "Capital Share over time", lw = 2, xlabel = L"$t$", ylabel = L"$\frac{r_t K_t}{Y_t}$", aspect_ratio = :none, ylimits=(0,1))
	p_kaldor4_ces = plot!(1:T₂-1, KS_ss_ces*ones(T₂-1), ls = :dash, lw=2, label = "Capital share (BGP)")

	plot(p_kaldor1_ces,p_kaldor2_ces, p_kaldor3_ces, p_kaldor4_ces, layout=(2,2))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
LaTeXStrings = "~1.3.0"
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "20323c8ab274ffda3db82f6f4180e906a358f762"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "9fb0b890adab1c0a4a475d4210d51f228bfc250d"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "7c29f0e8c575428bd84dc3c72ece5178caa67336"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.2+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "49cbf7c74fafaed4c529d47d48c8f7da6a19eb75"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.1"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "24b81b59bd35b3c42ab84fa589086e19be919916"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.11.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "47cf33e62e138b920039e8ff9f9841aafe1b733e"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.35.1+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═58753ecf-c322-4fb6-bb80-4abfc7cb5459
# ╟─1d082bd5-0086-4e29-aa5c-f4d41e1cf0fc
# ╟─0fee2ca0-7f57-4647-a187-36a9b3130cad
# ╟─f8f03434-6a0a-4c91-a084-7f86a13eaf93
# ╠═b2781e68-cc92-4186-a596-be753c52913e
# ╟─3ff3caf6-a1e3-4298-ac99-2ff430443dc7
# ╟─9e596114-e8ef-40ca-aeca-12d02bc49f21
# ╟─f70f5669-04bf-495a-8adc-8bc7afab6e2e
# ╟─66af3def-e62a-485a-ad93-901c3e920a96
# ╟─5669a9d4-c7cc-4426-a77d-33e23bc12175
# ╟─cc25bdab-c8df-4640-ae4f-d58ba44fc263
# ╟─29ee6c36-bcfb-40bc-bddf-7f8ba3d47c50
# ╟─be361a68-cf13-4493-a810-f327c5f97f2c
# ╠═3c9adefb-048c-4584-9ba8-b4ac22fdcd32
# ╟─53cc927f-02ec-4403-a3dc-6ea3791359b6
# ╟─66dbeb2a-32be-499a-a1d4-20c2647b8913
# ╠═2a5dd85f-8b53-4330-be69-3bee6b731caa
# ╟─76488fa6-2c94-4dbc-9c06-5bc4eaf6f3c6
# ╠═3e2db772-2806-4470-b308-b6a8cd271e46
# ╟─e417f1b3-fa3b-4736-80ef-d733fca46af8
# ╠═97f40165-4cc3-4705-b19b-a73b5f10ff6d
# ╟─c31be4ed-5d4d-4d9b-bd86-94c6020a3712
# ╠═031e0b88-732d-4ce2-9cd6-018faff82e84
# ╟─0d6af4bc-011f-457a-b4d0-60eda96aa460
# ╟─02d9cee7-9ca7-48cf-8663-f7fbc23ebc00
# ╟─17307ffa-4f1c-46e1-bb77-42c978464b91
# ╟─5d8e400e-9d1c-447c-8d4a-245d4f2a0fec
# ╟─b4c56d87-3b1c-4a06-974b-87f733bc8018
# ╟─04a3af31-73a1-4de0-9f65-6ab2e4721447
# ╟─9f6b6938-c148-48d3-ba78-39d6abd85403
# ╟─30332b4c-1efb-4236-8284-076a36686164
# ╟─e4e7a385-d0f8-4648-a097-0c3b74c46e0c
# ╟─56f6341b-47df-4c86-b503-2999feb34ddc
# ╟─900c31c6-12eb-4f4d-95c8-9aa4ec227beb
# ╟─f48174cb-706a-4966-bf8f-20d07f288712
# ╟─a7beeea9-44a9-432b-b41c-32c67ba92afd
# ╟─6932442c-a3dd-4dfb-b5ed-1bc14ec39c7b
# ╟─77287e52-6e9f-4228-8542-49fdffa9833f
# ╟─4ccdc633-e8eb-44bf-9d9f-95226b795ae3
# ╟─7fff17c8-b267-4d6f-b191-082ac3ac1a49
# ╟─a279af3a-db5c-4f64-9036-a3679328896c
# ╟─5ccb172d-c2ee-4b94-b9d3-3b44b059dd1a
# ╟─a486a61e-eb4d-4ad8-b7a6-fa3d95eeb682
# ╟─0203b0b7-958a-42e3-a064-d601f165b5fc
# ╟─643c3a17-5536-4ed6-8faa-f5980ebd70e6
# ╟─1ec69369-e4a8-4e51-9925-bf502204350f
# ╟─3cdd278a-1cb6-4ab7-92fb-7798a7eab02f
# ╟─9620ee95-69ca-44e3-892b-394b68a00a19
# ╟─00650bbb-f830-45e9-ab38-00c034ac8c2b
# ╟─ea573837-22e3-4f85-8ab6-08cd048dc1b1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
