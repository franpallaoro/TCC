getData = function(fdate, ldate){
  if (!require(BatchGetSymbols)) install.packages('BatchGetSymbols')
  tickers = c('OIBR4', 'TIMS3', 'VIVT3', 'GRND3', 'LAME4', 
              'BTOW3', 'COGN3', 'CVCB3', 'CYRE3', 'EZTC3', 'HGTX3', 
              'JHSF3', 'LCAM3', 'LREN3', 'MGLU3', 'MRVE3', 'RENT3', 
              'VVAR3', 'YDUQ3', 'ABEV3', 'BRFS3', 'BEEF3', 
              'JBSS3', 'MRFG3', 'PCAR3', 'CCRO3', 'EMBR3', 'GOLL4', 
              'POMO4', 'WEGE3', 'ECOR3', 'BBAS3', 'BBDC3', 'BBDC4', 
              'ITSA4', 'ITUB4', 'PSSA3', 'B3SA3', 'BBSE3', 'BRML3', 
              'CIEL3', 'IGTA3', 'MULT3', 'SANB11', 'SULA11', 'BRAP4', 
              'BRKM5', 'CSNA3', 'DTEX3', 'GGBR3', 'GGBR4', 'GOAU4', 
              'KLBN4', 'SUZB3', 'UNIP6', 'USIM5', 'VALE3', 'CESP6', 
              'CMIG3', 'CMIG4', 'CPLE6', 'EGIE3', 'ELET3', 'ELET6', 
              'LIGT3', 'SBSP3', 'TRPL4', 'CPFE3', 'ENBR3', 'ENEV3', 
              'ENGI11', 'EQTL3', 'TAEE11', 'PETR3', 'PETR4', 'UGPA3', 
              'CSAN3', 'PRIO3', 'FLRY3', 'HYPE3', 'QUAL3', 'RADL3')
  my_tickers <- paste0(tickers,'.SA')
  
  
  first_date <- fdate
  last_date <- ldate
  thresh_bad_data <- 0.99   
  bench_ticker <- '^BVSP'   
  
  
  
  l_out <- BatchGetSymbols(tickers = my_tickers,
                           first.date = first_date,
                           last.date = last_date,
                           bench.ticker = bench_ticker,
                           thresh.bad.data = thresh_bad_data)
  data = l_out$df.tickers[,7:9]
  data = BatchGetSymbols::reshape.wide(data)
  data = data$ret.adjusted.prices
  data = data[-which(is.na(data$ABEV3.SA)),]
  data = data[,-1]
  df = matrix(0, nrow(data), ncol(data))
  for (i in 1:ncol(data)) {
    df[,i] = data[,which(names(data) == my_tickers[i])]
  }
  df = as.data.frame(df)
  names(df) = my_tickers
  if(sum(is.na(df)) == 0){
    return(df)
  }else{
    return("NA nos dados")
  }
}

