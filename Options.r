CalcD1 <- function(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility) {
#
# Calculates the d1 input for BlackScholes model.
#
# Args:
#	See calcOptionPrice() for arg details
#
# Returns:
#	d1.

	d1 <- (log(stockPrice / strikePrice) + (costOfCarry + annVolatility ^ 2 / 2) *
		yearsToExpiration) / (annVolatility * sqrt(yearsToExpiration))
	
	d1
}


CalcD2 <- function(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility) {
#
# Calculates the d2 input for BlackScholes model.
#
# Args:
#	See calcOptionPrice() for arg details
#
# Returns:
#	d2.

	d2 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility) -
		annVolatility * sqrt(yearsToExpiration)
	
	d2
}


#Option price
CalcOptionPrice <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut = 'call') {
#
# Calculates the price of a call or put option according to the Black-Scholes model.
#
# Args:
#	callOrPut: a flag denoting whether it's a put or call. Only accepts 'call' or 'put'
#	stockPrice: the current stock price
#	strikePrice: the strike price
#	riskFreeRate: the risk free rate. Maturity should be equivalent to option's expiration.
#	yearsToExpiration: the number of years until option expiration. A decimal if within the year.
#	costOfCarry: the cost of carry. Assign this equal to the risk free rate if there are no
#				carry considerations. if a continuous dividend yield exists, this should be
#				equal to the risk free rate minus the dividend yield. If this is for a currency
#				option, then this should be equal to the risk free rate minus the risk free
#				rate times foreign risk free rate (r-rf).
#	annVolatility: the annualized volatility of the stock price
# 
# Returns:
#	The price of the option
	
	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	
	if (callOrPut == 'call') {
		optionPrice <- stockPrice * exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
			pnorm(tempD1, mean = 0, sd = 1) - strikePrice * exp(-riskFreeRate * yearsToExpiration) * 
			pnorm(tempD2, mean = 0, sd = 1)
	} else if (callOrPut == 'put') {
		optionPrice <- strikePrice * exp(-riskFreeRate * yearsToExpiration) * 
			pnorm(-tempD2, mean = 0, sd = 1) - stockPrice * exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
			pnorm(-tempD1, mean = 0, sd = 1)
	} else {
		optionPrice <- NULL
	}
	
	optionPrice
	
}
  

#Delta  
CalcDelta <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut = 'call') { 
#
# Calculates the delta of a call or put option.
# Delta is the first derivative of option value with respect to the price of the underlying.
# Expressed as the amount of money per underlying share that the option's value will gain or 
# lose as volatility rises or falls by 1%.
#
# Args:
#	See calcOptionPrice() for arg details
#
# Returns:
#	The Delta of the option  

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	
	if (callOrPut == 'call') {
		delta <- exp((costOfCarry - riskFreeRate) * yearsToExpiration) * pnorm(tempD1, mean = 0, sd = 1)
	} else if (callOrPut == 'put') {
		delta <- -exp((costOfCarry - riskFreeRate) * yearsToExpiration) * pnorm(-tempD1, mean = 0, sd = 1)
	} else {
		delta <- NULL
	}
	
	delta
}


#Vega
CalcVega <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) { 
#
# Calculates the vega of an option. 
# Vega is the first derivative of option value with respect to volatility of the underlying.
# Expressed as the amount of money per underlying share that the option's value will gain or 
# lose as volatility rises or falls by 1%.
#
# Args:
#	See calcOptionPrice() for arg details
#	Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Vega of the option  

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	
	vega <- stockPrice * exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
		dnorm(tempD1, mean = 0, sd = 1) * sqrt(yearsToExpiration)
		
	vega
}


#Theta
CalcTheta <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut = 'call') { 
#
# Calculates the Theta of a call or put option.
# Theta is the first derivative of option value with respect to the passage of time.
# This version is expressed as the amount of money per share of the underlying that the 
# option loses in one day (i.e. true theta / 365)
#
# Args:
#	See calcOptionPrice() for arg details
#
# Returns:
#	The Theta of the option  

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	
	if (callOrPut == 'call') {
		theta <- (-exp((costOfCarry - riskFreeRate) * yearsToExpiration) * stockPrice * 
			dnorm(tempD1, mean = 0, sd = 1) * annVolatility) / (2 * sqrt(yearsToExpiration)) -
			riskFreeRate * strikePrice * exp(-riskFreeRate * yearsToExpiration) * 
			pnorm(tempD2, mean = 0, sd = 1) + (riskFreeRate - costOfCarry) * stockPrice * 
			exp((costOfCarry - riskFreeRate) * yearsToExpiration) * pnorm(tempD1, mean = 0, sd = 1)
	} else if (callOrPut == 'put') {
		theta <- (-exp((costOfCarry - riskFreeRate) * yearsToExpiration) * stockPrice * 
			dnorm(tempD1, mean = 0, sd = 1) * annVolatility) / (2 * sqrt(yearsToExpiration)) + 
            riskFreeRate * strikePrice * exp(-riskFreeRate * yearsToExpiration) * 
			pnorm(-tempD2, mean = 0, sd = 1) - (riskFreeRate - costOfCarry) * stockPrice * 
			exp((costOfCarry - riskFreeRate) * yearsToExpiration) * pnorm(-tempD1, mean = 0, sd = 1)
	} else {
		theta <- NULL
	}
	
	theta / 365
}	


#Rho
CalcRho <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut = 'call') {
#
# Calculates the Rho of a call or put option.
# Rho is the first derivative of option value with respect to changes in the risk free rate.
# Expressed as the amount of money, per share of the underlying, that the value of the option 
# will gain or lose as the risk free rate rises or falls by 1.0% per annum (100 basis points).
#
# Args:
#	See calcOptionPrice() for arg details
#
# Returns:
#	The Rho of the option  

	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	if (callOrPut == 'call') {
		rho <- strikePrice * yearsToExpiration * exp(-riskFreeRate * yearsToExpiration) * 
			pnorm(tempD2, mean = 0, sd = 1)
	} else if (callOrPut == 'put') {
		rho <- -strikePrice * yearsToExpiration * exp(-riskFreeRate * yearsToExpiration) * 
			pnorm(-tempD2, mean = 0, sd = 1)
	} else {
		rho <- NULL
	}

	rho
}


#Gamma
CalcGamma <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Gamma of an option.
# Gamma is a second order derivative of option value with respect to changes in underlying price.
# Measures the rate of change of the option Delta with respect to changes in the underlying price.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Gamma of the option

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	optionGamma <- (exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
		dnorm(tempD1, mean = 0, sd = 1)) / (stockPrice * annVolatility * sqrt(yearsToExpiration))

	optionGamma
}


#Vanna	
CalcVanna <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Vanna of an option.
# Vanna is a second order derivative of option value, once to the underlying spot price and 
# once to volatility.
# AKA DdeltaDvol.
# Measures the sensitivity of the option Delta with respect to change in volatility.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Vanna of the option

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	vanna <- -exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
		dnorm(tempD1, mean = 0, sd = 1) * (tempD2 / annVolatility)

	vanna
}


#Charm
CalcCharm <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut = 'call') {
#
# Calculates the Charm of a call or put option.
# Charm is a second order derivative of option value, once to the underlying price and 
# once to time.
# AKA Delta decay, or DdeltaDtime.
# Measures the sensitivity of the option Delta with respect to passage of time.
#
# Args:
#	See calcOptionPrice() for arg details
#
# Returns:
#	The Charm of the option
                                    
	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	if (callOrPut == 'call') {
		charm <- -(costOfCarry - riskFreeRate) * exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
			pnorm(tempD1, mean = 0, sd = 1) + exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
			dnorm(tempD1, mean = 0, sd = 1) * (2 * costOfCarry * yearsToExpiration - tempD2 * 
			annVolatility * sqrt(yearsToExpiration)) / (2 * yearsToExpiration * annVolatility * 
			sqrt(yearsToExpiration))
	} else if (callOrPut == 'put') {
		charm <- (costOfCarry - riskFreeRate) * exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
			pnorm(-tempD1, mean = 0, sd = 1) + exp((costOfCarry - riskFreeRate) * yearsToExpiration) *
			dnorm(tempD1, mean = 0, sd = 1) * (2 * costOfCarry * yearsToExpiration - d2 * annVolatility * 
			sqrt(yearsToExpiration)) / (2 * yearsToExpiration * annVolatility * sqrt(yearsToExpiration))
	} else {
		charm <- NULL
	}

	charm
}


#Vomma
CalcVomma <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Vomma of an option.
# Vomma is a second order derivative of option value with respect to volatility
# AKA Vega Convexity.
# Measures the sensitivity of the option's Vega to changes in volatility
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Vomma of the option
               
	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	vomma = CalcVega(stockPrice, strikePrice, yearsToExpiration, riskFreeRate, costOfCarry, annVolatility) * 
		tempD1 * tempD2 / annVolatility

	vomma
}


#Veta
CalcVeta <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Veta of an option.
# Veta is a second order derivative of option value; once to volatility and once to time.
# AKA DvegaDtime.
# Measures the rate of change in Vega with respect to the passage of time.
# This version is divided by # of days per year to give the output the form: 
# change in Vega per one day.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Veta of the option

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	veta = CalcVega(stockPrice, strikePrice, yearsToExpiration, riskFreeRate, costOfCarry, annVolatility) * 
		((riskFreeRate - costOfCarry) + (costOfCarry * tempD1) / (annVolatility * sqrt(yearsToExpiration)) - 
		((1 + tempD1 * tempD2) / 2 * yearsToExpiration))
	
	veta / (365)	
}


#Speed
CalcSpeed <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Speed of an option.
# Speed is a third order derivative of option value wrt to underlying spot price.
# AKA Gamma of the Gamma.
# Measures the rate of change in Gamma with respect to changes in underlying price.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Speed of the option
	
	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	speed <- -(CalcGamma(stockPrice, strikePrice, yearsToExpiration, riskFreeRate, costOfCarry, 
		annVolatility) / stockPrice) * (tempD1 / (annVolatility * sqrt(yearsToExpiration)) + 1)
		
	speed
}


#Zomma
CalcZomma <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Zomma of an option.
# Zomma is a third order derivative of option value: twice to underlying price, once to volatility.
# AKA DgammaDvol.
# Measures the rate of change in Gamma with respect to changes in volatility.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Zomma of the option

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	zomma <- CalcGamma(stockPrice, strikePrice, yearsToExpiration, riskFreeRate, costOfCarry, 
		annVolatility) * ((tempD1 * tempD2 - 1) / annVolatility)
	
	zomma
}


#Colour
CalcColour <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Colour of an option.
# Colour is a third order derivative of option value: twice to underlying price, once to time.
# AKA Gamma decay.
# Measures the rate of change in Gamma with respect to the passage of time.
# This version is expressed as Gamma/year. Divide result by 365 to obtain Gamma/day.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Colour of the option

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	Colour <- -exp((costOfCarry - riskFreeRate) * yearsToExpiration) * 
		dnorm(tempD1, mean = 0, sd = 1) / (2 * stockPrice * yearsToExpiration * annVolatility * 
		sqrt(yearsToExpiration)) * (2 * (riskFreeRate - costOfCarry) * yearsToExpiration + 1 + 
		tempD1 * ((2 * costOfCarry * yearsToExpiration - tempD2 * annVolatility * 
		sqrt(yearsToExpiration)) / (annVolatility * sqrt(yearsToExpiration))))
	
	Colour
}


#Ultima
CalcUltima <- function(stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, annVolatility, callOrPut) {
#
# Calculates the Ultima of an option.
# Ultima is a third order derivative of option value wrt to volatility.
# AKA DvommaDvol.
# Measures the rate of change in Vomma with respect to changes in volatility.
#
# Args:
#	See calcOptionPrice() for arg details
#   Accepts callOrPut as a dummy variable strictly for consistency across the CalcGreek functions
#
# Returns:
#	The Ultima of the option

	tempD1 <- CalcD1(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)
	tempD2 <- CalcD2(stockPrice, strikePrice, yearsToExpiration, costOfCarry, annVolatility)

	ultima <- (-CalcVega(stockPrice, strikePrice, yearsToExpiration, riskFreeRate, costOfCarry, 
		annVolatility, callOrPut) / (annVolatility ^ 2)) * (tempD1 * tempD2 * 
		(1 - tempD1 * tempD2) + tempD1 ^ 2 + tempD2 ^ 2)
	
	ultima
}


#Volatility
CalcVolatility <- function(priceSeries, frequency = 'daily') {
#
# Calculates the annualized volatility of a stock's logarithmic returns.
#
# Args:
# 	priceSeries: a time series of stock prices
#	frequency: a flag for price frequency. May be one of the following:
#		daily
#		monthly
#		annual
#
# Returns:
#	A single value representing the annualized volatility.

	if (frequency == 'daily') {
		numPeriods <- 252
	} else if (frequency == 'monthly') {
		numPeriods <- 12
	} else if (frequency == 'annual') {
		numPeriods <- 1
	} else {
		stop("Incorrect frequency selected")
	}
	
	
	#transform the price series in to log returns	
	returnSeries <- diff(log(priceSeries))
	
	volatility <- sd(returnSeries, na.rm = TRUE)*sqrt(numPeriods)
	
	volatility	
}


#Implied Volatility
CalcImpVol <- function(optionPrice, stockPrice, strikePrice, yearsToExpiration, 
  riskFreeRate, costOfCarry, callOrPut = 'call'){
#
# Calculates the implied volatility of an option.
# This is an implementation of the bisection method. 
# See: http://en.wikipedia.org/wiki/Bisection_method
#
# Args: 
#	optionPrice: the market price of the option
#	See calcOptionPrice() for more arg details
#
# Returns:
#	The implied volatility of the option
	
	if (any(is.na(c(optionPrice, stockPrice, strikePrice, yearsToExpiration, 
		riskFreeRate, costOfCarry)))) {
		return(NA)
	}
	
	tolerance <- 0.0001
	lowerBoundVol <- 0.0001
	upperBoundVol <- 10
	impliedVol <- 0.25 #starting estimate	
	
	i=0
	while (abs(optionPrice - CalcOptionPrice(stockPrice, strikePrice, yearsToExpiration, 
		riskFreeRate, costOfCarry, impliedVol, callOrPut)) > tolerance) {
	
		i = i+1 
		if (i > 30) return(NA) #if this happens, the price is probably theoretically incorrect
		
		if (CalcOptionPrice(stockPrice, strikePrice, yearsToExpiration, 
		riskFreeRate, costOfCarry, impliedVol, callOrPut) < optionPrice) {
			lowerBoundVol <- impliedVol
		} else {
			upperBoundVol <- impliedVol
		}
		impliedVol <- mean(c(lowerBoundVol, upperBoundVol))
	}	
	impliedVol
}


# Option Expiration Date
ExtractOpExpiry <- function(optionCode) {
#
# Extracts the expiry date given an option code from getOptionChain using quantmod
#
# Args:
#	optionCode: a string containing a standardized option code. The format should
#				be similar to GOOG131129C00805000. i.e. Ticker, Yr, Mo, Da, C or P, 
#				Strike (00805), Decimal (000)
#
# Returns:
# A date object containing the expiration date of the option

	# Remove the ticker from the option code.
	digitsOfCode <- gsub("[^0-9]","",optionCode) 
	
	# Ignore the mini-option classifier of a '7' at the beginning of the digits
	ifelse(substr(digitsOfCode,1,1) == "7",
		dateString <- substr(digitsOfCode,2,7),
		dateString <- substr(digitsOfCode,1,6))
	
	expirationDate <- as.Date(dateString, "%y%m%d")
	
	expirationDate
}


# Option Data Frame
OptionDF <- function(optionChain, stockPrice, riskFreeRate, costOfCarry, strikeWindow = 0.25) {
#
# Transforms a list object returned by the quantmod getOptionChain() function by reducing
# the number of strikes to a more relevant window, and adding the implied volatilities
# for the bid and ask prices of each option
#
# Args:
#	optionChain: a list object returned from getOptionChain() function in quantmod package
#	riskFreeRate: the risk free rate. Currently a single number. Will eventually automate this
#	costOfCarry: see CalcOptionPrice for more details
#	strikeWindow: the +/- percentage range around the current stock price you want the 
#				  strike prices to be.
#
# Returns:
#	A list of the same format as the list in the getOptionChain() return. Contains 
#	strike prices, bids, asks, and implied volatilites for the bids and asks for 
#	all currently available options.


	ticker <- optionChain[[1]]$symbol
	optionDF <- data.frame()
	#Determine upper and lower bounds of the strike prices you want to look at
	strikeUB <- stockPrice * (1 + strikeWindow)
	strikeLB <- stockPrice * (1 - strikeWindow)
	
	#Loop through the expiraiton months (the 'names' of the chain)
	for (i in (1:length(names(optionChain)))) {
		#define the parts of the full option chain you'll be using
		calls <- optionChain[[i]]$calls[optionChain[[i]]$calls$Strike >= strikeLB & 
									     optionChain[[i]]$calls$Strike <= strikeUB,]
		puts <- optionChain[[i]]$puts[optionChain[[i]]$puts$Strike >= strikeLB & 
									   optionChain[[i]]$puts$Strike <= strikeUB,]
		slimChain = rbind(calls, puts)
		
		#pull relevant data from the option codes	
		optionCodeData <- SplitOptionCode(rownames(slimChain))
		
		#calculate the implied vols for bids and asks
		Bid.IVol <- mapply(CalcImpVol,slimChain$Bid, stockPrice, slimChain$Strike,
					      optionCodeData$YearsToExpiry, riskFreeRate, costOfCarry, 
					      optionCodeData$Type) 
		Ask.IVol <- mapply(CalcImpVol,slimChain$Ask, stockPrice, slimChain$Strike,
					      optionCodeData$YearsToExpiry, riskFreeRate, costOfCarry, 
					      optionCodeData$Type)
					
		#throw everything together
		combinedData <- cbind("optionCode" = rownames(slimChain),
							  optionCodeData,
							  "bid" = slimChain$Bid,
							  "ask" = slimChain$Ask, 
							  "bid.IVol" = round(Bid.IVol,6), 
							  "ask.IVol" = round(Ask.IVol,6))
		
		optionDF <- rbind(optionDF, combinedData)
	}
	
	#sort the output by expiry, then option type, then strike
	optionDF <- optionDF[order(optionDF$Expiry, optionDF$Type, optionDF$Strike),]
	rownames(optionDF) <- seq(1:nrow(optionDF))
	optionDF <- cbind(optionDF, stockPrice)
	optionDF
}


# Split apart option code
SplitOptionCode <- function(optionCode) {
#
# Picks apart an equity option code and returns its components.
#
# Args:
#	optionCode: A vector of equity option codes in string format.
#
# Returns:
#	A data frame with the following data: ticker, type of option, date of expiry,
#	years left until expiry, strike price, and a flag for mini-option status.

	#split apart the letters and numbers in the option code 
	codeLetters <- gsub("[^A-Z]","",optionCode)
	#If it's a mini-option, flag it, and adjust the numbers in the code for this edge-case
	codeNumbers <- gsub("[^0-9]","",optionCode)
	miniOptionFlag <- ifelse(substr(codeNumbers,1,1) == "7", "TRUE","FALSE")
	codeNumbers <- ifelse(substr(codeNumbers,1,1) == "7", substr(codeNumbers, 2, nchar(codeNumbers)),
					codeNumbers)
			
	callOrPut <- ifelse(substr(codeLetters, nchar(codeLetters), nchar(codeLetters)) == "C",
					"call", "put")
	
	ticker <- substr(codeLetters, 1, nchar(codeLetters)-1)
	
	expirationDate <- as.Date(substr(codeNumbers,1,6), "%y%m%d")
	
	yearsToExpiry <- as.integer(expirationDate - Sys.Date()) / 365
	
	strike <- as.integer(substr(codeNumbers,7,11)) + 
				as.integer(substr(codeNumbers,12,13)) / 100
		
	codeResults <- data.frame(
		ticker = ticker,
		type = callOrPut,
		miniOption = miniOptionFlag,
		expiry = expirationDate,
		yearsToExpiry = yearsToExpiry,
		strike = strike
		)
		
	codeResults
}


# Option payoff at expiration
OptionPayoff <- function(stockPrice, strikePrice, optionType) {
#
# Option payoff at expiration.
#
# Args:
#	optionType: Either 'call' or 'put'
#
# Returns:
#	The option payoff at expiration

	if(!(optionType == 'call' || optionType == 'put')) {
		return(NA)
	}
	
	ifelse(optionType == 'call', max((stockPrice - strikePrice), 0), max((strikePrice - stockPrice), 0))
}

# Make the output from quantmod's getOptionChain() more useable
UnlistChain <- function(optionList) {
#
# Converts the list object from getOptionChain into a simple data frame.
#
# Args:
#	optionList: output from getOptionChain() function
#
# Returns:
#	A data frame with option codes as row names, and the bid and ask prices

	chainDF <- unlist(optionList, recursive = FALSE, use.names = FALSE) #remove one list dimension (there are 2)
	chainDF <- chainDF[lapply(chainDF, is.character)==FALSE] #remove the elements of the list with ticker name
	chainDF <- do.call(rbind,chainDF) #combine all elements of list into dataframe
	chainDF <- chainDF[,c("Bid", "Ask")] #just choose the bid/ask prices
	chainDF
}	

