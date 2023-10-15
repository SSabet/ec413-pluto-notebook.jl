### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 58753ecf-c322-4fb6-bb80-4abfc7cb5459
begin
	#using Plots
	using PlotlyLight
	using LaTeXStrings
	#pyplot()
	Preset.Template.plotly_dark!()
	nothing
end

# ╔═╡ 1d082bd5-0086-4e29-aa5c-f4d41e1cf0fc
md"# _Solow Growth Model in Discrete Time (EC413-AT-2023)_

This is an **interactive** _Julia_ notebook to let you play around with examples of Solow model that you have seen in the lecture and seminars. 

**The good news is, you do not need to know _Julia_ to make use of this! For example, to see the effect of a change in parameters, change the parameter value and then press \"Shift+Enter\"**

Let's start by the Solow growth model defined in Assignment 2. A Key building block of the Solow model is the aggregate production function, so first we define that. You can change it to experiment with other functional forms if you wish (and press Shift+Enter for changes to take effect.)
"

# ╔═╡ 0fee2ca0-7f57-4647-a187-36a9b3130cad
md"""
## Convergence
### Assignment 2, Q1
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
#### Transforming the Model and Finding the Steady State of the Transformed Model
Letting $k_t$ denote capital per effective unit of labour and combining the transformed law of motion of capital with the production function we get

$k_{t+1}=\frac{1}{(1+g)(1+n)}\left[(1-\delta)k_t +s  (k_t)^\alpha \right].$

Steady state corresponds to $k_t$ not changing over time, i.e., $k_{t+1}=k_t = \bar{k}$:

$\bar{k}=\left(\frac{s}{g+n+gn+\delta}\right)^\frac{1}{1-\alpha}.$
"""

# ╔═╡ b9668385-f8b6-40a1-ba23-74f0e97dac92
md"""
#### Parameters
To get a better sense of dynamics of the model and be able to do plotting and numerical experimentation, we need to parameterise the model. Let's use the parameters from question 1 of assignment 2 (and the accompanying Excel file) for this purpose.

$\begin{align*}
n &= 0\\
g &= 0.2\\
\delta &= 0.5\\
\alpha &= 0.3\\
s &= 0.2
\end{align*}$
"""

# ╔═╡ adc5b813-3b34-4e38-994f-9a5a0d2449d4
# defining parameters; you can change them to see how other variables change
begin
	n = 0.0;
	g = 0.2;
	δ = 0.5;
	α = 0.3;
	s = 0.2;
	nothing
end

# ╔═╡ 8a3282c4-2353-43f0-8649-2f9abf4f126a
md"""
Now we can find the steady state value of $\bar{k}$ in the transformed model:
"""

# ╔═╡ bd8699e0-ebbe-4c61-83a4-c8291c1e262c
begin
	# finding the steady state value of ̄k in the transformed model
	k_ss = (s/(g+n+g*n+δ))^(1/(1-α))
	println(k_ss)
end

# ╔═╡ 7fff17c8-b267-4d6f-b191-082ac3ac1a49
md"""
#### Computing the Evolution of the Economy
Let's also choose the initial value of A, N and K:
"""

# ╔═╡ 90f25d99-5c59-4a2b-8f2c-d91955d59f41
begin
	A₁ = 1;
	N₁ = 1;
	K₁ = 0.1;
	nothing
end

# ╔═╡ 44135bf4-80d0-4b21-9669-0048095d944e
md"""
We can now compute the output in the first period using the aggregate productin function:
"""

# ╔═╡ d7b068b7-5ae7-4adb-8f5c-bfa52d7af6f5
Y₁ = F(A₁, K₁, N₁, α)

# ╔═╡ d40e886d-5a63-4547-a510-f7b602eb4436
md"""
Also choose the number of periods we want to simulate:
"""

# ╔═╡ a7acf601-6f39-4b4d-8890-ed98fde80e84
T = 10

# ╔═╡ 64192910-4c46-4238-bce1-c8bb6f6083f5
md"""
Now we start from the initial values of A, K and N and use the laws of motion to compute future values of these. Moreover we can use the aggregate production function to compute the output in each period.

Here is a function to do this. You can ignore it unless you want to peek into the sausage factory!
"""

# ╔═╡ df9e36a3-0ede-490b-bce2-9ad314092f98
function trajectories(A₁, K₁, N₁, T)
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

# ╔═╡ 93d415e3-10c4-4f70-8a6c-5e7903a11e09
md"""
Note that having all the aggrgate variables in hand, computing income per worker, capital per worker, capital per effective unit of labour etc. is straightforward:
"""

# ╔═╡ 569e746b-b7fc-439a-ab6e-ee94093e4bbc
begin
	A, K, N, Y = trajectories(A₁, K₁, N₁, T)
	YN = Y./N
	KN = K./N
	k = K./(A.*N)
end

# ╔═╡ 8ed00c7d-717b-4f9b-b61f-7697f0801970
md"""
#### Plotting the Evolution of England's Economy
Let's assume that England's starts from an initial level of capital per effective units of labour equal to half of Scotland (the steady state $\bar{k}$). We know that $K_t = A_t N_t \bar{k}_t$, so for England, $K_1$ is equal to:
"""

# ╔═╡ 5aa2bc72-34b2-409c-884c-2c41b6421271
K₁ₑ = A₁*N₁*(k_ss/2)

# ╔═╡ 13b4f021-01ca-41d4-992d-105e658fa79d
md"""
Now let's compute and plot trajectories of some variables of interest for England over time
"""

# ╔═╡ ae2fabc0-21dc-4b0f-bede-5e372a5c6635
begin
	Aₑ, Kₑ, Nₑ, Yₑ = trajectories(A₁, K₁ₑ, N₁, T)
	YNₑ = Yₑ./Nₑ
	KNₑ = Kₑ./Nₑ
	kₑ = Kₑ./(Aₑ.*Nₑ)
end

# ╔═╡ fe02f112-740c-4e98-9a04-42fcc2c0a7d4
md"""
Let's first have a look at the evolution of capital per effective worker $k$ for England (subscript e)
"""

# ╔═╡ 8732f2ed-2011-4d3e-8847-dccee208e10c
begin
	#plot(1:T, kₑ, label = "over time", lw = 2, xlabel = "time", ylabel = "kₜ")
	#plot!(1:T, k_ss*ones(T), label = "steady state", ls=:dash, lw = 2)
	p0 = Plot()(
	    x = 1:T, y = kₑ, name = "England"
	)(
	    x = 1:T, y = k_ss*ones(T), name = "Scotland (Steady State/BGP)"
	)
	
	p0.layout.title.text = "England: capital per effective unit of labour"
	p0.layout.xaxis.title.text = "time"
	p0.layout.yaxis.title.text = L"$k_t$"
	
	p0
end

# ╔═╡ ce64c5ba-9afe-4c8a-a854-8bfa0b25c09f
md"""
Let's also look, e.g., at how $\Delta k_{t+1} = k_{t+1} - k_t$ changes as function of $k_t$ as $k_t$ approaches the steady state
"""

# ╔═╡ 6516af4d-6a7f-46e1-ab9e-4721522ce79c
begin
	#	plot(kₑ[1:(T-1)], kₑ[2:T]-kₑ[1:(T-1)], label = "over time", lw = 2, xlabel = L"$k_t$", ylabel = L"$\Delta k_{t+1}$")
	p1 = Plot(x=kₑ[1:(T-1)], y= kₑ[2:T]-kₑ[1:(T-1)])
	p1.layout.title.text = "England"  # Make changes
	p1.layout.xaxis.title.text = L"$k_{t}$"
	p1.layout.yaxis.title.text = L"$\Delta k_{t+1}$"
	#p1.layout.xaxis_title.text = "xaxis"
	p1
end

# ╔═╡ 06ce56d1-6bdb-42ec-b302-b136dd1e3131
md"""
Now let's compute the same variables for Scotland (which starts from the steady state value of $k_1$) and plot the capital and income per capita of both countries
"""

# ╔═╡ 6fce058a-69c0-4f13-802a-2ca03346812d
begin
	K₁ₛ = A₁*N₁*k_ss
	
	Aₛ, Kₛ, Nₛ, Yₛ = trajectories(A₁, K₁ₛ, N₁, T)
	YNₛ = Yₛ./Nₛ
	KNₛ = Kₛ./Nₛ
	kₛ = Kₛ./(Aₛ.*Nₛ)
end

# ╔═╡ e24272cc-9d45-4679-9406-4fe3e64a67e6
begin
	#plot(1:T, Kₑ, label = "England", lw = 2, xlabel = "time", ylabel = L"$K_t$")
	#plot!(1:T, Kₛ, label = "Scotland", ls=:dash, lw = 2)
	
	p2 = Plot()(
	    x = 1:T, y = Kₑ, name = "England"
	)(
	    x = 1:T, y = Kₛ, name = "Scotland (Steady State/BGP)"
	)
	
	p2.layout.title.text = "England vs Scotland: capital (K)"  
	p2.layout.xaxis.title.text = "time"
	p2.layout.yaxis.title.text = L"$K_t$"
	
	p2
end

# ╔═╡ 56768d2a-6fa1-42e0-af6d-47e82fc379d1
begin
	#plot(1:T, YNₑ, label = "England", lw = 2, xlabel = "time", ylabel = L"$\frac{Y_t}{N_t}$")
	#plot!(1:T, YNₛ, label = "Scotland", ls=:dash, lw = 2)

	p3 = Plot()(
	    x = 1:T, y = YNₑ, name = "England"
	)(
	    x = 1:T, y = YNₛ, name = "Scotland (Steady State BGP)"
	)
	
	p3.layout.title.text = "England vs Scotland: income per capita (Y/N)" 
	p3.layout.xaxis.title.text = "time"
	p3.layout.yaxis.title.text = L"$\frac{Y_t}{N_t}$"
	
	p3
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlotlyLight = "ca7969ec-10b3-423e-8d99-40f33abb42bf"

[compat]
LaTeXStrings = "~1.3.0"
PlotlyLight = "~0.7.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "46e2d2a545189f983e7a431e8598cd71e2a6e72c"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Cobweb]]
deps = ["DefaultApplication", "Markdown", "OrderedCollections", "Random", "Scratch"]
git-tree-sha1 = "49e3de5be079f856697995001c587db8605506a9"
uuid = "ec354790-cf28-43e8-bb59-b484409b7bad"
version = "0.5.2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EasyConfig]]
deps = ["JSON3", "OrderedCollections", "StructTypes"]
git-tree-sha1 = "d22224e636afcb14de0cb5a0a7039095e2238aee"
uuid = "acab07b0-f158-46d4-8913-50acef6d41fe"
version = "0.1.15"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "95220473901735a0f4df9d1ca5b171b568b2daa3"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.13.2"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.PlotlyLight]]
deps = ["Artifacts", "Cobweb", "DefaultApplication", "Downloads", "EasyConfig", "JSON3", "Random", "Scratch", "StructTypes"]
git-tree-sha1 = "b842129ce0bc5fc230df8968738681a0db50223b"
uuid = "ca7969ec-10b3-423e-8d99-40f33abb42bf"
version = "0.7.4"

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

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"
"""

# ╔═╡ Cell order:
# ╟─58753ecf-c322-4fb6-bb80-4abfc7cb5459
# ╟─1d082bd5-0086-4e29-aa5c-f4d41e1cf0fc
# ╟─0fee2ca0-7f57-4647-a187-36a9b3130cad
# ╟─f8f03434-6a0a-4c91-a084-7f86a13eaf93
# ╠═b2781e68-cc92-4186-a596-be753c52913e
# ╟─3ff3caf6-a1e3-4298-ac99-2ff430443dc7
# ╟─9e596114-e8ef-40ca-aeca-12d02bc49f21
# ╟─b9668385-f8b6-40a1-ba23-74f0e97dac92
# ╠═adc5b813-3b34-4e38-994f-9a5a0d2449d4
# ╟─8a3282c4-2353-43f0-8649-2f9abf4f126a
# ╠═bd8699e0-ebbe-4c61-83a4-c8291c1e262c
# ╟─7fff17c8-b267-4d6f-b191-082ac3ac1a49
# ╠═90f25d99-5c59-4a2b-8f2c-d91955d59f41
# ╟─44135bf4-80d0-4b21-9669-0048095d944e
# ╠═d7b068b7-5ae7-4adb-8f5c-bfa52d7af6f5
# ╟─d40e886d-5a63-4547-a510-f7b602eb4436
# ╠═a7acf601-6f39-4b4d-8890-ed98fde80e84
# ╟─64192910-4c46-4238-bce1-c8bb6f6083f5
# ╟─df9e36a3-0ede-490b-bce2-9ad314092f98
# ╟─93d415e3-10c4-4f70-8a6c-5e7903a11e09
# ╠═569e746b-b7fc-439a-ab6e-ee94093e4bbc
# ╟─8ed00c7d-717b-4f9b-b61f-7697f0801970
# ╠═5aa2bc72-34b2-409c-884c-2c41b6421271
# ╟─13b4f021-01ca-41d4-992d-105e658fa79d
# ╠═ae2fabc0-21dc-4b0f-bede-5e372a5c6635
# ╟─fe02f112-740c-4e98-9a04-42fcc2c0a7d4
# ╟─8732f2ed-2011-4d3e-8847-dccee208e10c
# ╟─ce64c5ba-9afe-4c8a-a854-8bfa0b25c09f
# ╟─6516af4d-6a7f-46e1-ab9e-4721522ce79c
# ╟─06ce56d1-6bdb-42ec-b302-b136dd1e3131
# ╠═6fce058a-69c0-4f13-802a-2ca03346812d
# ╟─e24272cc-9d45-4679-9406-4fe3e64a67e6
# ╟─56768d2a-6fa1-42e0-af6d-47e82fc379d1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
