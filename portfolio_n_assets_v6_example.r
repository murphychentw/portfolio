library(PortfolioAnalytics)
library(quantmod)
library(PerformanceAnalytics)
library(zoo)
library(plotly)
library(ROI)
library(DEoptim)
library(doParallel)
require(ROI.plugin.glpk)
require(ROI.plugin.quadprog)

registerDoParallel()

sample_number <- 200

calculate_frontier <- function(returns.data) {

assets_number <- ncol(returns.data)
assets_name = colnames(returns.data)

returns.data <- na.omit(returns.data)

# Save mean return vector and sample covariance matrix
meanReturns <- colMeans(returns.data)
covMat <- cov(returns.data)

# Start with the names of the assets
port <- portfolio.spec(assets = assets_number, category_labels = assets_name)

# Box
port <- add.constraint(port, type = "box", min = 0.01, max = 1.0)
 
# Leverage
port <- add.constraint(portfolio = port, type = "weight_sum", min_sum=0.99, max_sum=1.01)

# Generate random portfolios
rportfolios_tmp <- random_portfolios(port, permutations = 5000, rp_method = "sample")
rportfolios = rportfolios_tmp[1:sample_number,]

# Get minimum variance portfolio
minvar.port <- add.objective(port, type = "risk", name = "var")

# Optimize
minvar.opt <- optimize.portfolio(returns.data, minvar.port, optimize_method = "random", rp = rportfolios)
								 
# Generate maximum return portfolio
maxret.port <- add.objective(port, type = "return", name = "mean")

# Optimize
maxret.opt <- optimize.portfolio(returns.data, maxret.port, optimize_method = "random", rp = rportfolios)

# Generate vector of returns
minret <- 0.06/100
maxret <- maxret.opt$weights %*% meanReturns
 
vec <- seq(minret, maxret[1], length = sample_number)

eff.frontier <- data.frame(Risk = rep(NA, length(vec)),
                           Return = rep(NA, length(vec)))
						   
frontier.weights <- mat.or.vec(nr = length(vec), nc = ncol(returns.data))
colnames(frontier.weights) <- colnames(returns.data)

for(i in 1:length(vec)){
  eff.port <- port
  eff.port <- add.constraint(eff.port, type = "return", name = "mean", return_target = vec[i])
  eff.port <- add.objective(eff.port, type = "risk", name = "var")
  eff.opt <- optimize.portfolio(returns.data, eff.port, optimize_method = "random", rp = rportfolios)
  
  eff.frontier$Risk[i] <- sqrt(t(eff.opt$weights) %*% covMat %*% eff.opt$weights)
  
  eff.frontier$Return[i] <- eff.opt$weights %*% meanReturns
  
  frontier.weights[i,] = eff.opt$weights
  
  print(paste(round(i/length(vec) * 100, 0), "% done..."))
}

cbind(eff.frontier, frontier.weights)
}

returns.data <- read.table("data_3_assets.txt")
returns.data$Annuity <- NULL
eff.frontier.A <- calculate_frontier(returns.data)

returns.data <- read.table("data_3_assets.txt")
eff.frontier.B <- calculate_frontier(returns.data)

p <- plot_ly(x = eff.frontier.A$Risk, y = eff.frontier.A$Return, 
        mode = "markers", type = "scattergl", name = "A: Bond+Stock",
        marker = list(size = 5, opacity = 1.0, color = "#000099" )) %>% 
  
  add_trace(x = eff.frontier.B$Risk, y = eff.frontier.B$Return, mode = "markers", 
            type = "scattergl", name = "B: Bond+Stock+Annuity",
            marker = list(color = "#990000", size = 5, opacity = 1.0)) %>% 
			
  layout(title = "Random Portfolios with Plotly",
         yaxis = list(title = "Mean Returns", tickformat = ".2%"),
         xaxis = list(title = "Standard Deviation", tickformat = ".2%"),
         plot_bgcolor = "#FFFFFF",
         paper_bgcolor = "#FFFFFF"
         )
p

p.weightA <- plot_ly() %>%
  add_trace(x = seq(1, sample_number), y = eff.frontier.A$Bond, 
        type = "bar", name = "Bond") %>% 
  
  add_trace(x = seq(1, sample_number), y = eff.frontier.A$Stock,
            type = "bar", name = "Stock") %>% 
			
  layout(title = "Portfolio Weighting", barmode = "stack",
         yaxis = list(title = "Weighting", tickformat = ".2%"),
         xaxis = list(title = "Index"),
         plot_bgcolor = "#FFFFFF",
         paper_bgcolor = "#FFFFFF"
         )
p.weightA
		 
p.weightB <- plot_ly() %>%
  add_trace(x = seq(1, sample_number), y = eff.frontier.B$Annuity,
            type = "bar", name = "Annuity") %>% 

  add_trace(x = seq(1, sample_number), y = eff.frontier.B$Bond, 
        type = "bar", name = "Bond") %>% 

  add_trace(x = seq(1, sample_number), y = eff.frontier.B$Stock,
            type = "bar", name = "Stock") %>% 	

  layout(title = "Portfolio Weighting", barmode = "stack",
         yaxis = list(title = "Weighting", tickformat = ".2%"),
         xaxis = list(title = "Index"),
         plot_bgcolor = "#FFFFFF",
         paper_bgcolor = "#FFFFFF"
         )
p.weightB
