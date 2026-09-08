library(testthat)
library(dssatutils)
fixture <- testthat::test_path('..','fixtures','agera5_2001_complete.csv')

test_that('year bounds require every boundary day, including leap day', {
  p <- tempfile(fileext='.WTH'); on.exit(unlink(p))
  dates <- seq(as.Date('2000-01-01'),as.Date('2000-12-31'),by='day')
  rows <- paste(format(dates,'%Y%j'),'12 25 15 2 10 60 3')
  writeLines(rows,p); expect_true(is_wth_valid(p,2000,start_year=2000))
  for (part in list(rows[-1],head(rows,-1),rows[183],rows[-60])) {
    writeLines(part,p); expect_false(is_wth_valid(p,2000,start_year=2000))
  }
  writeLines(rows[1:60],p)
  expect_true(is_wth_valid(p,2000,start_year=2000,end_date='2000-02-29'))
  expect_false(is_wth_valid(p,2000,start_year=2000,end_date='2000-03-01'))
})

test_that('cache validates individual forcing values and multicell calendars', {
  p<-tempfile(fileext='.csv'); on.exit(unlink(p))
  original<-read.csv(fixture,check.names=FALSE)
  for (bad in c(NA,Inf,-99)) {
    x<-original; x$Solar_Radiation_Flux[6]<-bad; write.csv(x,p,row.names=FALSE)
    expect_false(dssatutils:::.agera5_validate_timeseries_csv(p,2001))
  }
  other<-original;other$longitude<- -89.9
  write.csv(rbind(original,other),p,row.names=FALSE)
  expect_true(dssatutils:::.agera5_validate_timeseries_csv(p,2001))
  expect_false(dssatutils:::.agera5_validate_timeseries_csv(p,2001,area=c(40.05,-90.05,39.95,-89.95)))
  write.csv(rbind(original,other[-1,]),p,row.names=FALSE)
  expect_false(dssatutils:::.agera5_validate_timeseries_csv(p,2001))
})

test_that('R download accepts staged response despite post-transfer exception', {
  skip_if_not_installed('ecmwfr')
  work<-tempfile();dir.create(work);on.exit(unlink(work,recursive=TRUE))
  local_mocked_bindings(.agera5_ensure_ecmwfr_key=function(...) invisible(TRUE), .package='dssatutils')
  local_mocked_bindings(wf_request=function(request,path,...) {
    file.copy(fixture,file.path(path,request$target));stop('post-transfer error')
  },.package='ecmwfr')
  job<-list(year=2001,area=c(40.05,-90.05,39.95,-89.95),cache_dir=work,data_format='csv')
  dest<-dssatutils:::.agera5_download_timeseries_job(job)
  expect_true(file.exists(dest))
  expect_true(dssatutils:::.agera5_validate_timeseries_csv(dest,2001))
  expect_equal(dssatutils:::.agera5_download_timeseries_job(job),dest)
})

test_that('parallel R assembly reads complete cached years without credentials or network', {
  work<-tempfile();dir.create(work);on.exit(unlink(work,recursive=TRUE))
  local_mocked_bindings(.agera5_ensure_ecmwfr_key=function(...) stop('unexpected credentials'), .package='dssatutils')
  points<-data.frame(ID='00000001',LAT=40,LONG=-90)
  area<-dssatutils:::.agera5_split_timeseries_chunks(40,-90,.1)[[1]]$area
  for(y in 2001:2002) {
    x<-read.csv(fixture,check.names=FALSE); x$valid_time<-sub('2001',as.character(y),x$valid_time)
    write.csv(x,dssatutils:::.agera5_timeseries_cache_path(work,y,area),row.names=FALSE)
  }
  out<-file.path(work,'weather')
  process_weather_agera5(points,2001,2002,out,'ID','LAT','LONG',2,NULL,work,
                        agera5_backend='timeseries',agera5_timeseries_chunk_degrees=.1,cache_only=TRUE)
  expect_true(is_wth_valid(file.path(out,'00000001.WTH'),2002,start_year=2001))
})
