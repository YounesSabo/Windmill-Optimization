
	using Markdown
	using InteractiveUtils
	using Pkg
		cd(dirname(dirname(@__DIR__)))
		Pkg.activate(pwd())

	using Plots
	using CSV
	using DataFrames
	using DelimitedFiles
	using Optim
	using Logging
	using StatsPlots
	using Dates
	using LaTeXStrings
	using LoggingExtras
	

	# Kijken wat het aantal threads zijn voor het te laten runnen 
	Threads.nthreads()
	# we werken nu met 5 THREADS 

	# OPSTELLEN VAN DE LOGGER
	
		io = open(joinpath(@__DIR__,"LoggerStroming.txt"),"w+")
		logger = ConsoleLogger(io, Logging.Debug)
		global_logger(logger)
	

md"""
##  Vaste Constantes doorheen de simulatie
Binnen de simulatie zal er met een constante Ct en β gewerkt worden. We werken in een vak van 2000m op 2000m waarbij er een disretisatie van dx gelijk aan 0.5m is. 
- Ct = 0.75
- β = 0.1
- ρ = 1.225 [{kg}{m^3}]
De eigenschappen van de turbine is de volgende:
https://www.vestas.com/en/products/offshore/V236-15MW/V236-15MW
- Type = V236-15MW
- Diameter = 232 m
- S = 42.273,27 m^2

Verder zal men 4 verschillende turbines gebruiken waarbij enkel θ de veradnerlijke is. De posities worden ook vast gekozen op basis van de opstelling die in de realiteit wordt gebruikt.
- Turbine 1 : (500, 500) 
- Turbine 2 : (1500, 500) 
- Turbine 3 : (1500, 1500) 
- Turbine 4 : (500, 1500) 
"""


# vast waardes van het veld
	const X = 2000; #[m]
	const Y = 2000; #[m]
	const dx = 0.5;
	const dy = 0.5;

	const lenX = trunc(Int, X/dx)
	const lenY = trunc(Int, Y/dy)

# vaste waardes van de stroming
	Ct = 0.75
	β = 0.1
	ρ = 1.225
# vaste waardes van de turbines
	const Diameter = 232 #[m]
	const Diameter_disc = trunc(Int64,Diameter/dx)
	const S = π*(Diameter/2)^2
	const Cp = 1 # Alle wind door de turbine zal naar het vermogen gaan 
	# posities van de turbines
	# turbine 1
	PosT1 = (1000, 1000)
	# turbine 2
	PosT2 = (3000, 1000)
	# turbine 3
	PosT3 = (3000,3000)
	# turbine 4
	PosT4 = (1000, 3000)



md"""
## Data van de windsnelheden
Binnen deze simulatie zal de data van het KMNI voor de noordzee gebruikt worden tijdens 2020.

DATA eigenschappen

- HH = de tijd per uur
- DD = de windrichting
- FH = gemiddelde snelheid per uur
- FF = windsnelheid van de laatste 10min van het uur
- FX = Hoogste windstoot (in 0.1 m/s) over het afgelopen uurvak 
- T = De temperatuur
"""



	input = joinpath(@__DIR__,"gegwind_EURO_2020.csv")
	data = CSV.read(input,DataFrame)



md"""
### Functies om de data om te zetten
We zullen het gemiddelde per dag gaan nemen gezien we data per uur hebben voor elke dag. Daarnaast zullen we de windrichtingen indelen volgens het compas.

Omzetten van de windrichtingen 

	WindConv(W::Vector{Any})
Omzetten naar daggemiddeldes

	DayConv(FH)
"""


		function DayConv(FH)
			gem = 0
			arr = []
			for (i,v) in enumerate(FH)
					
				if mod(i,7*24) == 1
					push!(arr,gem/(7*24));
					gem = 0;
				end
				gem += v;
			end
			return arr[2:end]
		end


	# STRUCTUREN

	mutable struct Flow 
		u::Array{Float64,2}
		v::Array{Float64,2}
		α::Array{Float64,2}
		U0::Float64

		function Flow(U0)
			u = zeros(lenY,lenX)
			v = zeros(lenY,lenX)
			α = zeros(lenY,lenX)

			new(u,v,α,U0)
		end

	end


	struct Turbine 
			D::Int64
			θ#::Float64
			Pos::Tuple{Int64, Int64}
			function Turbine(D::Int, θ,Pos::Tuple{Int64,Int64}) #::Float64 weg bij θ
				
				new(D,θ,Pos)
			end 
	end
		Base.show(io::IO, t::Turbine) = println(io,"\nTurbine @$(t.Pos) with\n-diameter $(t.D)\n-inclination $(t.θ*180/π)°")


	struct TurbineField 
		
		field :: Array{Int,2}
		Turbines::Array
		WR::Int
		LocTurb :: Dict
		Δθ::Int

		function TurbineField(θ,wr,Pos = [PosT1,PosT2,PosT3,PosT4]) 
			
			if (typeof(θ) <: Vector{Float64}) == false 
				throw(error("θ must be <: Vector{Float64}, not $(typeof(θ)) when calling TurbineField"))
			elseif ((typeof(Pos)) <: Vector{Tuple{Int64, Int64}}) == false
				throw(error("Pos must be <: Vector{Tuple{Int64, Int64}}, not $(typeof(Pos)) when calling TurbineField"))
			end

			field = zeros(lenY,lenX)
			
			# Aantal turbines
			lenTurb = length(Pos)
			Δθ = 0
			LocTurb = Dict()
			

			# afhankelijk van waar de wind komt moeten we onze array gaan aanpassen 

			# EXTREMITEITEN NOORD EN ZUID
			if wr == 0
				# noord
				for t in 1:lenTurb
					y,x = Pos[t]
					
					for i in 0:trunc(Int64,Diameter_disc/2)
						field[y-i,x] = t
						field[y+i,x] = t
					end
				end

				field = transpose(field)
				
				for i in 1:lenTurb
					LocTurb[i] = Tuple.(findall(v->v==i,field))
				end
			elseif wr == 180
				# zuid
				for t in 1:lenTurb
					y,x = Pos[t]
			
					for i in 0:trunc(Int64,Diameter_disc/2)
						field[y-i,x] = t
						field[y+i,x] = t
						
					end
				end
				
				field = rotr90(field)
				
				for i in 1:lenTurb
					LocTurb[i] = Tuple.(findall(v->v==i,field))
				end
			end

			# hier willen we operaties doen voor w's die niet gelijk zijn aan 0 of 180
			if wr > 180
				# turbine 1 en 2 zullen de eerste turbines zijn die de wind zal zien
				if wr >315
					# de diagonaal, we roteren het veld over een hoek van -45 graden, deze hoek moet bij θ opgeteld worden
					Δθ = -45
					for t in 1:lenTurb
						y,x = Pos[t]
						
						for i in 0:trunc(Int64,Diameter_disc/2)
							field[y-i,x-i] = t
							field[y+i,x+i] = t
						end
					end

				elseif wr<225
					# de ANTIdiagonaal, we roteren het veld over een hoek van +45 graden, deze hoek moet bij θ opgeteld worden
					Δθ = 45
					for t in 1:lenTurb
						y,x = Pos[t]
						
						for i in 0:trunc(Int64,Diameter_disc/2)
							field[y-i,x+i] = t
							field[y+i,x-i] = t
						end
					end
				else
					for t in 1:lenTurb
						y,x = Pos[t]
				
						for i in 0:trunc(Int64,Diameter_disc/2)
							field[y-i,x] = t
							field[y+i,x] = t
							
						end
					end
				end

				for i in 1:lenTurb
					LocTurb[i] = Tuple.(findall(v->v==i,field))
				end

			elseif wr<180 && wr>0
				if wr <45
					# Eerst leggen we de turbine op de ANTIdiagonaal en daarna spiegelen we de matrix ==> komen uiteindeljk op de diagonaal
					Δθ = -45
					for t in 1:lenTurb
						y,x = Pos[t]
						
						for i in 0:trunc(Int64,Diameter_disc/2)
							field[y-i,x+i] = t
							field[y+i,x-i] = t
						end
					end

				elseif wr > 135
					# Eerst leggen we de turbine op de diagonaal en daarna spiegelen we de matrix ==> komen uiteindeljk op de ANTIdiagonaal
					Δθ = 45
					for t in 1:lenTurb
						y,x = Pos[t]
						
						for i in 0:trunc(Int64,Diameter_disc/2)
							field[y-i,x-i] = t
							field[y+i,x+i] = t
						end
					end
				else
					for t in 1:lenTurb
						y,x = Pos[t]
				
						for i in 0:trunc(Int64,Diameter_disc/2)
							field[y-i,x] = t
							field[y+i,x] = t
							
						end
					end
				end
				reverse!(field, dims = 2)
				for i in 1:lenTurb
					LocTurb[i] = Tuple.(findall(v->v==i,field))
				end
			end

			
			
			Turbines = [Turbine(Diameter,θ[i],Pos[i]) for i in 1:lenTurb]
			new(field,Turbines,wr,LocTurb,Δθ)

		end
	end
		Base.show(io::IO, t::TurbineField) = println(io,"Windmill parc with $(length(t.Turbines)) windmills.\nWind comming from $(t.WR)\n Turbines:\n$(t.Turbines)\n With field inclination $(t.Δθ)")



			# test om struct uit te proberen
				#Hoeken = [9.0,-9.0,8.0,9.0]*π/180
				#TurbField = TurbineField([9.0],270,[PosT1])


	# FUNCTIES VOOR HET BEREKENEN VAN DE STROMINGEN

	function ControlFlow(U::Array{Float64,2},V::Array{Float64,2},α::Array,IntTurb::Int64,U_init::Number) 
			
			u = deepcopy(U)
			A = deepcopy(α)

			# Inclinatie van de turbine
			θ = TurbField.Turbines[IntTurb].θ + TurbField.Δθ

			# De posities van de wieken van de turbine
			LocTurb = TurbField.LocTurb[IntTurb]

			# Windrichting van de wind
			wr = TurbField.WR

			# alles tot de turbine is U0 + herdefinieren van stroming.U0 voor het correcte model
			for i in 1:length(LocTurb)
				y = LocTurb[i][1]
				for x in 1:LocTurb[i][2]
					
					if wr>270
						diff = wr-270
						u[y,x] = abs(U_init*sin((TurbField.WR-diff)*π/180))
					elseif wr>180 && wr<270
						diff = wr-180
						u[y,x] = abs(U_init*cos((TurbField.WR-diff)*π/180))
					elseif wr>90 && wr<180
						diff = wr-180
						u[y,x] = abs(U_init*cos((TurbField.WR-diff)*π/180))
					else
						u[y,x] = U_init*sin((TurbField.WR)*π/180)
					end
					Stroming.U0 = u[y,x]
				end

				
			end

			
			U0 = Stroming.U0
			# vanop de turbine verandert het
			for (y,x) in LocTurb
				
				α_n = (cos(θ))^2*sin(θ)*Ct/2
				ΔU = (cos(θ))^3*U0*Ct/2
				U_n = U0-ΔU
				
				u[y,x] = U_n*cos(α_n)
				A[y,x] = α_n
			end

			# vanaf de turbine zal er een uitbreidende cylinder ontstaan met een andere snelheid

			# x zal beginnen vanaf de x van de turbine
			# Gezien een discretisatie hier gelijk is aan 0.2m zal er om de 10 vakjes een uitbreiding van de cylinder zijn. 

			Lim = Dict(1 => [-1], length(LocTurb) => [-1])
			
			Threads.@threads for i in 1:length(LocTurb)  
				y = LocTurb[i][1]
				
				for x in LocTurb[i][2]+1:size(u,2)
					
					# afstand tov turbine vinden 
					X_int = LocTurb[i][2]
					Dist = (x-X_int)*0.2

					# formules toepassen
					δ = Diameter+β*Dist
					
					α_n = θ == 0 ? 0 : (Diameter/δ)^2*(cos(θ))^2*sin(θ)*Ct/2
					ΔU = (Diameter/δ)^2*(cos(θ))^3*U0*Ct/2
					U_n = U0-ΔU
					
					# inbrengen van het schuin effect
					# Kiezen 1/α_n aangezien deze waarde steeds groter zal worden naarmate we verder van de turbine gaan, zo kunnen we het exponentieel convergereend effect vormen.
					
					if 	α_n != 0 && mod(x-X_int, trunc(Int,0.5/α_n)) == 0 && y>1 && y < 3999
						y += trunc(Int,1*sign(θ)) 
						
						u[y,x] = U_n*cos(α_n)
						A[y,x] = α_n
						
					end	
									
					#UITBREIDING
					Expand = 10
					

					if i == 1
						
						if mod(x-X_int,Expand) < 10^-5 && abs(Lim[i][end]) < y-1 
							newval = Lim[i][end] - 1
							push!(Lim[i],newval)
						end

						for dir in Lim[i]
							if y+dir>0
							u[y+dir,x] = U_n*cos(α_n)
							A[y+dir,x] = α_n
							end
						end
						
					elseif i == length(LocTurb)
						
						if mod(x-X_int,Expand) < 10^-5 && y + Lim[i][end] < size(u,1)-1
							newval = Lim[i][end] + 1
							push!(Lim[i],newval)
						end

						for dir in Lim[i]
							# Ook hier weer checken of de uitbreiding onze array verlaat.
							if dir + y < size(u,1)
							u[y+dir,x] = U_n*cos(α_n)
							A[y+dir,x] = α_n
							end
						end
						
					end

					u[y,x] = U_n*cos(α_n)
					A[y,x] = α_n

				end
				
			end
			
			
			
			
			
		return u,A
	end


	function CalcFlow(flows,stroming,PosTurb)
		#FUNCTIE DIE STROMING.U ZAL OSPTELLEN
		# Moeten alle stromen uitrekenen van de turbines die voor de laatste turbines staan
		First = first.(PosTurb)
		n = count(v->v == First[1],First) != length(PosTurb) ? filter!(e->e ∉ findall(e->e==findmax(First)[1],First),[i for i in eachindex(PosTurb)] ) : [i for i in eachindex(PosTurb)] # het moet een vector worden!
		
	
		for i in n
			Un,αn = ControlFlow(stroming.u,stroming.v,stroming.α,i,stroming.U0)
			flows[i] = (Un,αn)
		end

		# We gaan het gehele veld opstellen voor de stormingen
		

		for i in n
			stroming.u += flows[i][1]
		end

		for y in 1:size(stroming.u,1),x in 1:size(stroming.u,2)
			if stroming.u[y,x] == 0
				stroming.u[y,x] = stroming.U0
			end
		end	
		
	end


	# BEREKNEN VAN HET VERMOGEN
	# P = \frac{1}{2} \cdot ρ \cdot Cp \cdot S \cdot V^3 Met Cp = 1

	function CalcPower(Loc::Dict,U::Array{Float64,2})
		Power = zeros(length(Loc))
		
		for i in 1:length(Loc)
			LocTurb = Loc[i]

			N = length(LocTurb)
			
			Utot = 0
			for (y,x) in LocTurb
				Utot += U[y,x-1]
			end
			
			Ugem = Utot/N
			
			P = 0.5*ρ*S*Ugem^3
			Power[i] = P
		end
		
		return sum(Power)
	end


	# BEREKENEN VAN HET VERMOGEN IN EEN BEPAALD WEEKINTERVAL MET GEGEVEN θ

	function PowerWeek(θ,optim = false,weeks::Array = Weken,PosTurb::Vector{Tuple{Int64, Int64}} = [PosT1,PosT2,PosT3,PosT4])
	
		Power = zeros(length(weeks))
		@warn("Calculation for the power are starting @time $(now())")
		
		for w in eachindex(weeks) #Weken
			i = weeks[w]
			U0 = Vmean[i]
			
			# TurbField moet over heel het document gekend zijn. 
			global TurbField = TurbineField(θ, Windrichting_num[i],PosTurb)
			Flows = Dict()
			global Stroming = Flow(U0)

			# Uitrekenen van het veld
			CalcFlow(Flows,Stroming,PosTurb)

			# Uitrekenen van het vermogen
			Power[w] = CalcPower(TurbField.LocTurb , Stroming.u)
		end
		
		# gezien we met optim.jl werken moeten we de negatieve output minimaliseren, dit komt overeen met het maximaliseren van de positieve output
		if optim
			
			return -1*sum(Power) # gebruiken we voor de optimalisatie
		end
		
		@info("Total power = $(sum(Power))")
		flush(io)
		return Power #abs.(Power)  
	end



# DATA VOOR DE BEREKENINGEN
 const Vmean = DayConv(data.:("FH,"))
 const Windrichting_num = round.(Int,DayConv(data.:("DD,")))	
 const Weken = [i for i in 1:length(Vmean)]

# BESPREKING MODEL 
# - Ideale situatie van de luchtstroom ==> gaan turbulentie vermijden dus geen overlapping van stromen
# - Turbulentie zal niet gesimuleerd worden ==> geen LES of RANS
# - Actuator disk model wordt gebruikt
# - Cp = 1 ==> Alle wind zal gebruikt worden om naar vermogen omgezet te worden
# - Enkel rotatie van de gehele turbine wordt geoptimaliseerd


# KLEINSCHALIGE DEMO
	# Voorbeeld waar we 1 turbine gaan optimaliseren
	# We stellen de turbine op PosT1 op voor een θ en een windrichting
	wr = 270 #[°]
	θ_test = [3*π/180] #[rad]
	TestTurb = TurbineField(θ_test,wr,[PosT1])

	# Berkenen de jaarlijkse output voor deze ene turbine
	Output_Test = PowerWeek(θ_test,false,Weken,[PosT1])

	function VisTest(res)
		# Zullen de weekelijkse output voorstellen
		# Zullen de procentuele stijging in power weergeven 
		p1 = Plots.bar(res,widen = false,color=:blue,label=:"Weekly Power",legend=:topright );

		increase = (cumsum(res) - pushfirst!([i for i in cumsum(res[1:end-1])],0))/sum(res)	
		p2 = Plots.bar(increase);
		savefig(p1,joinpath(@__DIR__,"Grafieken/Test/WeaklyOutput.png"))
		savefig(p2,joinpath(@__DIR__,"Grafieken/Test/Increment.png"))
	end
	
	VisTest(Output_Test)
	
	# Optimaliseren van het veld met 1 turbine voor een maximale output tijdens de winter
		lx = ones(length(θ_test))*-10*π/180
		ux = ones(length(θ_test))*10*π/180
		WinterWeeks = [i for i in 1:11];push!(WinterWeeks,51);push!(WinterWeeks,52)

		obj_functionWINTER(x) = PowerWeek(x,true,WinterWeeks,[PosT1])

		options = Optim.Options(store_trace = true,show_trace = true, iterations = 5, outer_iterations = 3) 
		@info("Optimization for the winter weeks is starting")
		res_winter = Optim.optimize(obj_functionWINTER,lx,ux,θ_test,Fminbox(NelderMead()),options)

		MINIMIZER_test = res_winter.minimizer*180/π #*180/π voor omzetten naar degrees
		MINIMUM = res_winter.minimum

	# jaarlijkse output van de optimalisatie
		Output_OptimTest = PowerWeek(MINIMIZER_test,false,Weken,[PosT1])
		VisTest(Output_OptimTest)
		sum(Output_OptimTest)
	# We krijgen een optimale hoek dat net geen 0° is, dit klopt gezien men de kleinste afnamen van de snelheid heeft bij 0° graden.





# ANALYSE VAN DE RESULTATEN
	function Optimfun(fun::Function,Weeks::Vector,inner_it::Int64,outer_it::Int64,θ_init)
		@warn("Optimization for the winter weeks is starting @time $(now())")
		obj_functionWINTER(x) = fun(x,true,Weeks)

		# Indien het aantal itteraties niet gespecificeerd is zal er automatisch 1000 itteraties worden uitgevoerd
		lx = ones(length(θ_init))*-10*π/180
		ux = ones(length(θ_init))*10*π/180
		options = Optim.Options(store_trace = true,show_trace = true, iterations = inner_it, outer_iterations = outer_it) 
		println("Optimization for the winter weeks is starting with:\n-outer itterations = $(outer_it)\n-inner itterations = $(inner_it)\n-for function $(fun)")

		t1 = convert(Dates.Millisecond,now())
		res_winter = Optim.optimize(obj_functionWINTER,lx,ux,θ_init,Fminbox(NelderMead()),options)
		t2 = convert(Dates.Millisecond,now())

		ComputationTime = t2-t1
		@info ComputationTime

		MINIMIZER = res_winter.minimizer*180/π #*180/π voor omzetten naar degrees
		
		@info res_winter
		@info(" *MINIMIZER\n    $MINIMIZER") 

		# Jaarlijkse output berekenen van de optimale hoeken
		OUTPUT_year = PowerWeek(MINIMIZER,false,Weken)

		flush(io)
		return MINIMIZER,OUTPUT_year,res_winter
		
	end
	θ_init = [1,0,1,0]*π/180
	# WINTER OPTIMALISATIE
		
	WinterWeeks = [i for i in 1:11];push!(WinterWeeks,51);push!(WinterWeeks,52)
	MINIMIZER_winter,OUTPUT_winter,res_winter = Optimfun(PowerWeek,WinterWeeks,3,5,θ_init)

	# Het resultaat
	θ_winter = [9.823980322679526, -0.6864592832862105, -0.4278097863361681, 4.213306691416667]*π/180 #[°] in de array
	Output_optW = PowerWeek(θ_winter,false)
	sum(Output_optW)


	# JAARLIJKSE OPTIMALISATIE

	MINIMIZER_year,OUTPUT_year,res_year = Optimfun(PowerWeek,Weken,5,20,θ_init)
	θ_jaar =  [9.841288690931625, 2.7777230977179803, 0.5354271486581488, 9.480643325371258]
	Output_optY = PowerWeek(θ_jaar,false)
	sum(Output_optY)
	# Werkelijke cijfers (zal een maandelijkse weergave zijn) https://www.belgianoffshoreplatform.be/nl/productie-data/
	Output_real = [628.387,803.233,704.035,340.262,419.229,361.094,370.658,356.351,436.545,877.656,639.344,792.631]
	Output_realM = [628.387, 628.387, 628.387, 628.387, 803.233, 803.233, 803.233, 803.233, 704.035, 704.035, 704.035, 704.035, 704.035, 340.262, 340.262, 340.262, 340.262, 419.229, 419.229, 419.229, 419.229, 361.094, 361.094, 361.094, 361.094, 361.094, 370.658, 370.658, 370.658, 370.658, 356.351, 356.351, 356.351, 356.351, 436.545, 436.545, 436.545, 436.545, 436.545, 877.656, 877.656, 877.656, 877.656, 639.344, 639.344, 639.344, 639.344, 792.631, 792.631, 792.631, 792.631, 792.631]
	plot(Output_real/sum(Output_real))
# VISUALISATIE VAN DE RESULTATEN

	# WINTER
	let
		pw = Plots.bar(Output_optW,widen = false,ylabel=L"Vermogen [W]",xlabel=L"Week",
						legend=:topright,label="Opgewekt Vermogen",title="Wekelijkse vermogen voor een optimale winter",margin=7Plots.mm,grid = false)

		pw_inc = Plots.bar(Output_optW/sum(Output_optW),widen = false,ylabel=L"Stijging [\%]",
							xlabel=L"Week",legend=:topright,label="Wekelijkse stijging",title="Wekelijkse bijdrage tot het totaal vermogen",margin=7Plots.mm,grid=false)

		pw_tot = Plots.bar(cumsum(Output_optW),widen = false,ylabel=L"Vermogen [W]",xlabel=L"Week",legend=:topright,label="Cummulatief vermogen",title="Totaal Vermogen",margin=7Plots.mm,grid=false)
		Winter = Output_optW[1:11];push!(Winter,Output_optW[51]);push!(Winter,Output_optW[52])
		Lente = Output_optW[12:24]
		Zomer = Output_optW[25:37]
		Herfst = Output_optW[38:50]
		tot = sum(Output_optW)
		pw_seizoen = groupedbar([sum(Winter) sum(Lente) sum(Zomer) sum(Herfst)],
								bar_position = :stack,
								bar_width=0.7,
								label=["Winter $(round(sum(Winter)/tot,digits = 4)*100)%" "Lente $(round(sum(Lente)/tot,digits = 4)*100)%" "Zomer $(round(sum(Zomer)/tot,digits = 4)*100)%" "Herfst $(round(sum(Herfst)/tot*100,digits = 2))%"],
								orientation="h",
								yaxis = false,
								margin=7Plots.mm,
								title="Vermogen per seizoen",
								xlabel=L"Vermogen [W]")

		L = @layout [ [a{0.5w} b{0.5w}]
					c{0.3h}]
		subplots = plot(pw,pw_tot,pw_seizoen,; layout=L)
		# final (global plot)
		global_title = plot(title = "Optimalisatie voor de wintermaanden met een totaal vermogen van 1.4625e13 W ", grid=false, showaxis=false, ticks=false, bottom_margin = 10Plots.mm,upper_margin = 10Plots.mm)
		pw_final = plot(global_title, subplots, layout=@layout([A{0.01h}; B]) ,size = (1100,900))

		savefig(pw,joinpath(@__DIR__,"Grafieken/Winter/WeeklyPower.png"))
		savefig(pw_inc,joinpath(@__DIR__,"Grafieken/Winter/WinterIncrement.png"))
		savefig(pw_tot,joinpath(@__DIR__,"Grafieken/Winter/WinterTotaal.png"))
		savefig(pw_seizoen,joinpath(@__DIR__,"Grafieken/Winter/WinterSeizoen.png"))
		savefig(pw_final,joinpath(@__DIR__,"Grafieken/Winter/final.png"))
	end



	# JAAR
	let
		py = Plots.bar(Output_optY,widen = false,ylabel=L"Vermogen [W]",xlabel=L"Week",
						legend=:topright,label="Opgewekt Vermogen",title="Wekelijkse vermogen",margin=7Plots.mm,grid = false)

		py_inc = Plots.plot(Output_optY/sum(Output_optY),widen = false,ylabel=L"Stijging [\%]",
							xlabel=L"Week",legend=:topright,label="Wekelijkse stijging",title="Wekelijkse bijdrage tot het totaal vermogen",margin=7Plots.mm,grid=false);plot!(Output_realM/sum(Output_real),
							label = "Werkelijke maandelijkse bijdrage")
		py_tot = Plots.bar(cumsum(Output_optY),widen = false,ylabel=L"Vermogen [W]",xlabel=L"Week",legend=:topright,label="Cummulatief vermogen",title="Totaal Vermogen",margin=7Plots.mm,grid=false)
		WinterY = Output_optY[1:11];push!(WinterY,Output_optY[51]);push!(WinterY,Output_optY[52])
		LenteY = Output_optY[12:24]
		ZomerY = Output_optY[25:37]
		HerfstY = Output_optY[38:50]
		totY = sum(Output_optY)
		py_seizoen = groupedbar([sum(WinterY) sum(LenteY) sum(ZomerY) sum(HerfstY) ],
								bar_position = :stack,
								bar_width=0.7,
								label=["Winter $(round(sum(WinterY)/totY,digits = 4)*100)%" "Lente $(round(sum(LenteY)/totY,digits = 4)*100)%" "Zomer $(round(sum(ZomerY)/totY*100,digits = 2))%" "Herfst $(round(sum(HerfstY)/totY*100,digits = 2))%"],
								orientation="h",
								yaxis = false,
								margin=7Plots.mm,
								title="Vermogen per seizoen",
								xlabel=L"Vermogen [W]")

		L = @layout [ [a{0.5w} b{0.5w}]
					c{0.3h}]
		subplots = plot(py,py_tot,py_seizoen,; layout=L)
		# final (global plot)
		sum(WinterY)/totY
		global_title = plot(title = "Optimalisatie voor het hele jaar met een totaal vermogen van 1.62e13 W ", grid=false, showaxis=false, ticks=false, bottom_margin = 10Plots.mm,upper_margin = 10Plots.mm)
		py_final = plot(global_title, subplots, layout=@layout([A{0.01h}; B]) ,size = (1100,900))

		savefig(py,joinpath(@__DIR__,"Grafieken/Jaar/JaarPower.png"))
		savefig(py_inc,joinpath(@__DIR__,"Grafieken/Jaar/JaarIncrement.png"))
		savefig(py_tot,joinpath(@__DIR__,"Grafieken/Jaar/JaarTotaal.png"))
		savefig(py_seizoen,joinpath(@__DIR__,"Grafieken/Jaar/JaarSeizoen.png"))
		savefig(py_final,joinpath(@__DIR__,"Grafieken/Jaar/final.png"))
	end

	
# VERDERE BESPREKING
# -Een betere aanpak zou zijn om naar een CFD toe te gaan waar alle soorten effecten van de storming in betrekkingen genomen kunnen worden
