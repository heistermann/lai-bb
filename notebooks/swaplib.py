import pandas as pd
import numpy as np
import os
import re
import subprocess
import shutil
from io import StringIO
import requests, io
import zipfile
import scipy
import matplotlib.pyplot as plt

#PySWAP dependencies
from pathlib import Path

def q_as_mm(q, A, interval="Y"):
    """Discharge from m3/s to mm/a
    """
    if interval=="Y":
        x = 365.25 * 86400
    elif interval=="M":
        if q.index.freqstr==interval:
            x = q.index.days_in_month * 86400
        else:
            print("WARNING: data is not in monthly interval!")
            x = 30.5 * 86400
    else:
        print("Cannot deal with interval %s, yet." % interval)
        raise
    return 1000 * q * x / (A * 1e6)


def set_par_in_string(string, parname, value, printedits=False):
    """Replace a value of parameter parname in string by value.
    """
    #import re
    #string="HLIM3H = -600.0            ! Pressure head below which water uptake reduction starts at high Tpot [-1d4..100 cm, R]"
    #parname = "HLIM3H_F"
    is_factor = re.search(pattern='_f$', string=parname, flags=re.I) is not None
    if is_factor:
        parname = parname[:-2] #remove "_f" from parameter name
    #    # pattern = re.compile(r'^\s*%s\s*=' % parname)
    #     pattern = "^\s*%s\s*=\s*([\-\.0-9]+).*" % parname
    #     #pattern = "^\s*%s\s*=\s*([\-\.0-9]+).*" % "parname"
    #     old_value = re.sub(string=string, pattern=pattern, repl="\\1") #extract old parameter value from file
    #     try:
    #         old_value=float(old_value)
    #         value = old_value * value   #the given value is interpreted as a factor to be combined with the original one found in the file
    #     except:
    #         print("WARNING: Old value for parameter %s could not be read, ignored." % parname)
    #         return string

    #items = re.findall("^.*%s*=*.*$" % parname,string,re.MULTILINE)
    #items = re.findall("^\s+%s\s+=*" % parname,string,re.MULTILINE)
    pattern = re.compile(r'^\s*%s\s*=' % parname)
    items = [line for line in string.split('\n') if pattern.search(line)]
    #print(items)
    for item in items:
        if item.strip()[0]=="*":
            items.remove(item)
    if len(items)==0:
        print("WARNING: no parameter %s found in string." % parname)
        return string
    if len(items)>1:
        print("AMBIGUITY ERROR: %d occurences of parameter %s:" % (len(items),parname))
        print(items)        
        sys.exit()
    item = items[0]
    #print(items)
    parval, comment = item.strip().split("!")
    par, old_value = parval.split("=")
    if is_factor:
        try:
            old_value = float(old_value)
            value = old_value * value  # the given value is interpreted as a factor to be combined with the original one found in the file
        except:
            print("WARNING: Old value for parameter %s could not be read, ignored." % parname)
            return string

    replacement = "  " + par.strip() + " = " + str(value)# + "   !" + comment
    replacement += (29-len(replacement))*" " + "!" + comment
    if printedits:
        print(replacement)
    return string.replace(item, replacement)


def replace_between(string, newtab, marker_start, marker_end="* End of table"):
    """Replace in string everything between marker_start and marker_end by newtab.
    """
    to_replace = string[string.find(marker_start):string.rfind(marker_end)]
    return string.replace(to_replace, newtab.to_string(index=False)+"\n")

def replace_table(string, newtab_arg, marker_end="* End of table"):
    """Replace a table in string by newtab.
    """

    if (newtab_arg.__class__ == pd.core.frame.DataFrame):
        key="_" #set _ as indicator that the key should not be used
        newtab = newtab_arg.copy() #use as is
        is_factor = np.any([re.search(pattern='_f$', string=s, flags=re.I) is not None for s in newtab_arg.columns]) #does ths table contain multipliers?
        if is_factor: #remove "_f" from column names
            newtab.columns = [re.sub(pattern='_f$', string=s, repl="", flags=re.I) for s in newtab_arg.columns]

    if (newtab_arg.__class__ == dict):
        key   = list(newtab_arg.keys())[0]
        newtab= list(newtab_arg.values())[0]
        is_factor = re.search(pattern='_f$', string=key, flags=re.I) is not None  # does this table contain multipliers?
    if is_factor and len(newtab.index) != 1:
        print("WARNING: Table %r containing multipliers can only have a single row." % newtab)
        return string

    if key[0] == "_":
        #pattern = "\s+".join([s + "(_f)*" for s in newtab.columns]) #use table columns for searching position
        pattern = "\s+".join([s + "" for s in newtab.columns]) #use table columns for searching position
    else:
        pattern = "\s*"+key+"\s*=\s*" #use variable name for searching position


    pattern = "("+pattern + ".*?)"+ re.escape(marker_end)
    #pattern = "("+pattern + ")(.*?)"+ re.escape(marker_end) #find also (variable) markerstring
    pattern = re.compile(pattern=pattern, flags=re.DOTALL+re.I)
    to_replace = pattern.findall(string)

    #with open('test.txt', 'w') as f:
    #    f.write(string)

    #print(to_replace)
    if len(to_replace)==0:
        print("WARNING: table %r not found in string." % newtab)
        return string
    if len(to_replace)>1:
        print("WARNING: Multiple replacements of table %r ." % newtab)

    if is_factor:
        #read old values from file
        tablestr = StringIO(to_replace[0])
        old_vals = pd.read_csv(filepath_or_buffer=tablestr, delim_whitespace=True) #read old values
        #newtab = newtab_arg
        #newtab.columns = ['CROPSTART', 'CROPEND', 'CROPFIL', 'CROPTYPE_f']
        #newtab = newtab.iloc[0:1]
        #multiply values
        for i in range(0, len(old_vals.columns)):
            pass
            fac = newtab.iloc[0,i]
            if ~ np.isnan(fac):
                old_vals.iloc[:,i] =  fac * old_vals.iloc[:,i][:]

        newtab = old_vals #set updated table

    #use space as column separator
    table_str = newtab.to_string(index=False, header= key[0] == "_")
    if key[0] != "_":
        table_str = "\n" + key + " =\n" +table_str   #prepend variable identifier
    string = string.replace(to_replace[0], table_str+"\n")
    #string = string.replace(to_replace[0][0], table_str+"\n")

    # preserve found marker string: problematic, as the separator in the following table needs to be the same
    # table_str = newtab.to_string(index=False, header= False)
    # table_str = re.sub(pattern=" ", repl="\t", string=table_str) #replace space with tab
    # #re.sub(pattern=pattern, repl="\\1 "+table_str, string=string)
    # string = re.sub(pattern=pattern, repl="\\1\n"+table_str+"\n"+marker_end, string=string)

    return string

def set_pars_in_file(infile, outfile, pardict, tables):
    """
    Sequence to replace control parameters in a SWAP input file.
    
    First replace parameter values, then tables.

    Parameters
    ----------
    infile: str
        name of input file used as template to replace the requested parameters, e.g. *.swp, *.crp, *.bbc
    outfile: str
        name of output file containing the modified parameters
    pardict: dict
        dictionary with keys identifying the parameters to be replaced in infile. Parameter names ending with "_f" will be interpreted as factors affecting the original values found in the file.
    tables: list of Pandas-dataframes,
        ...with correspondingly-named columns

    Returns
    -------
    nothing, write modified file to outfile

    Notes
    -----
    Tables to be replaced can be passed via argument "pardict" or argument "tables".
    When passed via "tables"  their position is identified via their column names, e.g. for ["ISUBLAY", "ISOILLAY", ...] for *.swp; ["DATE1", "GWLEVEL"] *.bbc

    When passed via "pardict", their position is identified their key, e.g. "dateharvest" for *.crp.
    When the key starts with _, the key is ignored and the position is identified via column names as if supplied via argument "tables".
    """


    with open(infile) as f:
        string = f.read()
    for key, val in pardict.items():
        if val.__class__ == pd.core.frame.DataFrame: #special case: replace a table
            string = replace_table(string, {key:val})
        else:
            string = set_par_in_string(string, key, pardict[key])
    for table in tables:
        string = replace_table(string, table)
    with open(outfile, "w") as f:
        f.write(string)


def modify_metfile(infile, outfile, add_cols=dict(), mult_cols=dict()):
    """
    Modify SWAP meteo-file input file additively or multiplicatively with specified offset or factor

    Parameters
    ----------
    infile: str
        input file to be modified, e.g. a *.met-file
    outfile: str
        output file containing the modified meteo-input for SWAP
    add_cols: dict, default dict()
        dictionary with keys identifying the variables (i.e. columns in met-file) to be altered by offset (additively)
    mult_cols: dict, default dict()
        dictionary with keys identifying the variables (i.e. columns in met-file) to be altered by offset (additively)

    Returns
    -------
    nothing, write modified file to outfile

    Notes
    -----
    """
    met_data = pd.read_csv(filepath_or_buffer=infile, sep=",") #read original met-data

    # alter column by given factor
    for key,val in mult_cols.items():
        if key not in met_data.columns:
            print("Could not find column "+key+" in meteo-file, ignored.")
            continue
        met_data[key] = met_data[key] * val

    # alter column by given offset
    for key,val in add_cols.items():
        if key not in met_data.columns:
            print("Could not find column "+key+" in meteo-file, ignored.")
            continue
        met_data[key] = met_data[key] + val

    met_data.to_csv(path_or_buf=outfile, sep=',', index=False)

def get_vgm_table(parlist):
    """Return a dataframe from list of VGM parameters
    
    Example: parlist=[0.05, 0.38, 0.0050, 1.550, 7.04, 0.5, 0.0, 7.04, 1300.0]
    """
    vgmtab = pd.DataFrame(
        columns=["ORES", "OSAT", "ALFA", "NPAR", "KSATFIT", "LEXP", "H_ENPR", "KSATEXM", "BDENS"],
        data = [parlist])
    return vgmtab


def get_vdiscr_table(maxdepth, depths=[0, 5, 15, 50, 100, 200, 500, 8000], 
                     stepsize=[  1.,2.5, 5., 10., 10., 20.,  50.]):
    """Create vertical discretization table for swp file.
    
    maxdepth: float, maximum depth of soil column (cm)
    depths : list of len(stepsize)+1, depth thresholds between which stepsize is applied
    stepsize: list of len(depths)-1, defines the discretization step between depths[i] and depths[i+1]
    """
    discr = pd.DataFrame(columns=["ISUBLAY", "ISOILLAY", "HSUBLAY", "HCOMP", "NCOMP"])
    #maxdepth = 1500
    #depths =   [0, 5, 15, 50, 100, 200, 500, 4000]
    #stepsize = [  1.,2.5, 5., 10., 10., 20.,  50.]
    ISUBLAY = 1
    for i in range(len(stepsize)):
        bottom = depths[i+1]
        top = depths[i]
        if top >= maxdepth:
            continue
        if bottom > maxdepth:
            bottom = maxdepth
        discr.loc[i] = [i+1, 1, bottom-top, stepsize[i], (bottom-top)/stepsize[i]]
    discr["ISUBLAY"] = discr["ISUBLAY"].astype(int)
    discr["ISOILLAY"] = discr["ISOILLAY"].astype(int)
    discr["NCOMP"] = discr["NCOMP"].astype(int)
    return discr

def get_crop_rotation_table(ystart, yend, cropfil, croptype):
    """Return DataFrame with SWAP crop rotation table
    """
    df = pd.DataFrame(columns=["CROPSTART", "CROPEND", "CROPFIL", "CROPTYPE"])
    df["CROPSTART"] = ["%d-01-01" % i for i in range(ystart, yend+1)]
    df["CROPEND"] = ["%d-12-31" % i for i in range(ystart, yend+1)]
    df["CROPFIL"]= "'%s'" % cropfil
    df["CROPTYPE"]= croptype
    return df


def read_swap_vap(vapfile="data/swap/result.vap", depths=None, var="wcontent"):
    """Read one variable var out of vapfile at depths.
    """
    with open(vapfile, encoding = "ISO-8859-1") as f:
        fout = f.read()
        fout = fout.replace(" ","")
    sim = pd.read_csv(StringIO(fout), skiprows=11)
    if not var in sim.columns:
        print("ERROR: var %s in not in %s." % (var, vapfile))
    sim = sim[["date", "top",var]]
    sim.columns = ["date", "depth",var]
    sim["depth"] = -sim.depth
    sim["date"] = pd.to_datetime(sim.date)
    #sim = sim.set_index("date")
    if depths is None:
        depths = np.unique(sim.depth)
    else:
        # check if target depths are actually present in dataset
        depths_in_data = np.in1d(depths, np.unique(sim.depth))
        if not np.all(depths_in_data):
            print("Not all depths are in %s: %r" % (vapfile, depths[depths_in_data]))
            sys.exit()
    depths = np.array(depths)#.astype("int")
    dates = np.unique(sim.date)
    dates.sort()
    sim2 = pd.DataFrame(index=dates, columns=depths)
    for depth in depths:
        try:
            tmp = sim.loc[sim.depth==depth]
            tmp = tmp.sort_values(by="date")
            sim2[depth] = tmp[var].to_numpy()
        except:
            print("Problem with depth %d" % depth)
            continue
    return sim2


def read_swap_default(f):
    """Read SWAP output file at path f and return DataFrame.
    """
    with open(f, encoding = "ISO-8859-1") as f:
        fout = f.read()
        fout = fout.replace(" ","")
    sim = pd.read_csv(StringIO(fout), comment="*")
    cols = np.array(sim.columns)
    # Unify date column label
    cols[cols=="date"] = "Date"
    sim.columns = cols.tolist()
    sim["Date"] = pd.to_datetime(sim.Date)
    sim = sim.set_index("Date")
    return sim

def run_swap(runid,
             savefiles=["data/swap/result.inc"],
             saveto="data/ezg/Nuthe/results/",
             silent=True, swap_path=None, rundir="./", rm_rundir=False):
    """Runs SWAP and saves "savefiles" under "saveto" using "runid" to name files.
    """
    if swap_path is None:
        swap_path = "./run_swap.sh" if not is_windows() else "./run_swap.bat"

    # run model
    if silent:
        stdout = stderr =  subprocess.DEVNULL
    else:
        stdout = stderr =  None

    # spout = subprocess.run([swap_path],
    #                        stdout = stdout,
    #                        stderr = stdout,
    #                        input="\n",   #allows termination of process, if something went wrong and SWAP gets stuck waiting for ENTER
    #                        cwd=rundir)
    result = _run_exe(tempdir=rundir, swap_path = swap_path) #call executable

    if 'normal completion' in result:
        try:
            for f in savefiles:
                file_name, file_extension = os.path.splitext(f)
                shutil.move(f, os.path.join(saveto, "%s%s" % (runid, file_extension)))
            result=0
        except:
            pass
            #result=-2
    else:
        #print(result)
        pass
        #result=-1
        # raise Exception(
        #     f'Model run failed. \n {result}')

    if rm_rundir and result == 0:
        shutil.rmtree(rundir, ignore_errors=True)  # delete rundir after completion, unless there had been an error
    return result

#copied from PySWAP for potential later use
def _run_exe(tempdir: Path, swap_path=None) -> str:
    if swap_path is None:
        swap_path = Path(tempdir, 'swap.exe') if is_windows() else './swap420'

    p = subprocess.Popen(swap_path,
                         stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         cwd=tempdir)

    return p.communicate(input=b'\n')[0].decode()

def is_windows() -> bool:
    from platform import system
    """Checks if the current OS is Windows."""
    return True if system() == 'Windows' else False


def list_files_in_url_directory(url):
    """List files in url directory.
    """
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.split("\n")
        files = []
        for line in lines:
            if "href=" in line:
#                print(line)
                start_index = line.find("href=") + 6
                end_index = line.find('"', start_index)
                files.append(line[start_index:end_index])
        return files
    else:
        print("Failed to fetch directory listing:", response.status_code)
        return []

# Global variable with columns in "produkt" files of DWD daily collectives
product_columns = {
    # climate
    "KL" : ["id", "datetime", "QN_3", "maxwindspeed", "avgwindspeed", "QN_4", "precip", "preciptype",
           "sunhours", "snowheight", "cloudiness", "vaporpress",  "press", "tempmean", "relhum", "tempmax",
           "tempmin", "tempmin5cm", "eor"],
    # precipitation
    "RR" : ["id","datetime","QN_6","precip","precipform","snowdepth","newsnowdepth","eor"],
    # radiation
    "ST" : ["id","datetime","QN_592","atmlong","solardiff","glorad","sunhours","eor"]
}

def get_station_from_collective(stationid, directory_url, collective, returnvals=["data"]):
    """Retrieves data from collective for stationid from DWD's opendata as DataFrame.

    Parameters
    ----------
    stationid: integer
        DWD station ID
    directory_url: str
        source directory on DWD-server, e.g."https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/more_precip/historical/"
    collective: str
        e.g. "RR" or "KL"
    returnval: list of str, default: ["data"]
        list, indicating if "data", "meta" or both are to be returned.

    Returns
    -------
    (tuple of) dataframes, containing "data", "meta" or both

    Notes
    -----
    """
    import re

    files_in_url = list_files_in_url_directory(directory_url)
    tmp = "%05d" % stationid
    df=None #for collecting the actual data
    meta = dict() #for collecting meta info
    for file in files_in_url:
        if not ("tageswerte_%s_" % collective.upper() + tmp) in file:
            continue
        url = directory_url+file
        r = requests.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        for inzip in z.namelist():
            #extract data
            if ("produkt" in inzip) and ("data" in returnvals):
                z.extract(inzip)
                df = pd.read_csv(inzip, sep=";", na_values=[-999,-9999])
                df.columns = product_columns[collective]
                df["datetime"] = pd.to_datetime(df.datetime, format="%Y%m%d")
                df = df.set_index("datetime")

            #extract metadata
            if (re.match(pattern="Metadaten.*\.txt$", string=inzip) is not None) and ("meta" in returnvals):
                z.extract(inzip)
                meta_var = re.sub(pattern=".*_([^_]*)_0*" + str(stationid) + ".*\\.txt$", repl="\\1", string=inzip) #extract name of meta-variable
                meta[meta_var]  = pd.read_csv(inzip, sep=";", na_values=[-999,-9999], skipfooter=1, engine="python", encoding='ansi') #get metadata

            try: #delete unzipped file, if existing
                os.remove(inzip)
            except:
                pass

        ret_val=list() #assemble return value
        if df is None:
            if ("data" in returnvals):
                print("No product found in %s" % file)
        else:
            ret_val = ret_val + [df]

        if len(meta)==0:
            if "meta" in returnvals:
                print("No metadata found in %s" % file)
        else:
            ret_val = ret_val + [meta]

        return(ret_val)



def plot_trend(starttrend, endtrend, df, varname, ax=None, **kwargs):
    """Plot trend line from starttrend to endtrend for column varname in DataFrame df.
    """
    tmptrend = df[starttrend:endtrend]
    tmptrend = tmptrend[~tmptrend[varname].isna()]
    tmptrend["delta"] = [(item - tmptrend.index[0]).total_seconds() for item in tmptrend.index]
    slope, intercept, r, p, se = scipy.stats.linregress(tmptrend.delta.to_numpy(), tmptrend[varname].to_numpy())
    #slope2 = q_as_mm(slope*365.25*86400,metadata.at[i,"Ae"])*10
    ypred = intercept + slope*tmptrend.delta
    if not ax is None:
        plt.sca(ax)
    plt.plot(tmptrend.index, ypred, **kwargs)


def read_metf(f):
    """Read SWAP met file to DatFrame.
    """
    metdf = pd.read_csv(f)
    metdf["datetime"] = ["%s-%s-%s" % (metdf.at[i,"YYYY"], metdf.at[i,"MM"], metdf.at[i,"DD"]) for i in metdf.index]
    metdf["datetime"] = pd.to_datetime(metdf.datetime)
    return metdf.set_index("datetime")


def resample_with_nan(df, freq, maxnan, fun="mean"):
    """Account for a maximum number of nans (maxnan) when resampling.
    """
    res = df.resample(freq).agg(fun)
    #notna = df.isna().resample(freq).sum()
    notna = df.isna().resample(freq).sum()
    res[notna>maxnan] = np.nan
    return res

def read_dwdvar(collective, varname, datadir="data/dwd/"):
    """Read variable varname out of DWD collective from local directory datadir.
    
    collective : string, one out of KL, RR, or ST
    varname : string, climate variable name contained in collective
    datadir : string, local path
    """
    f = os.path.join(datadir, collective, "%s-%s.csv" % (collective, varname))
    df = pd.read_csv(f)
    df["datetime"] = pd.to_datetime(df.datetime)
    df = df.set_index("datetime")
    df.columns = [int(i) for i in df.columns]
    return df

def predict_glorad(df, model, featnames):
    try:
        df = df.drop(columns=["ypred"])
    except:
        pass
    df1 = df[featnames]
    df2 = df1[~np.any(df1.isna(), axis=1)]
    df2["ypred"] = model.predict(df2.to_numpy())
    df = pd.concat([df,df2[["ypred"]]], axis=1)
    return df

def rmse(x, y):
    return np.sqrt(np.mean((x-y)**2))

def mae(x, y):
    return np.mean(np.abs(x-y))

def r2(x, y):
    return np.corrcoef(x, y)[0,1]**2


def retentioncurve(x, theta_R, theta_S, alpha, n, saturation_index=False):
    """Returns theta from h.
    """
    m = 1 - 1/n
    sat_index = (1 + (alpha * x)**n)**(-m)
    if (saturation_index):
        out = sat_index
    else:
        out = theta_R + (theta_S - theta_R) * sat_index
    return(out)

def Kr_theta_curve(theta, thetaS, thetaR, n, Ks, f=0.5):
    m = 1 - (1/n)
    Se = (theta - thetaR)/(thetaS - thetaR)
    a = (1 - (1 - Se**(1/m))**m)**2
    out = Ks * (Se**f) * a
    return(out)


def Kr_h_curve(h, alpha, n, Ks, f=0.5):
    m = 1 - (1/n)
    Se = (1/(1 + (alpha * h)**n))**m
    b = (1 - (1 - Se**(n/(n - 1)))**m)**2
    out = Ks * (Se**f) * b
    return(out)
    
def drop_gappy_gauges(df, start, end, minperc):
    """Removes columns of time series for which a minimum data availability is not given.
    
    Parameters
    ----------
    df : dataframe with datetime index
    start : datetime object or string to mark start of evaluation period
    end : datetime object or string to mark end of evaluation period
    minperc : minimum data availability (in percent) required from a column to be kept
    
    Returns
    -------
    dataframe : from which gappy columns are removed
    """
    dtimes = pd.date_range(start, end, freq="D")
    df = df.loc[start:end]
    if len(dtimes) < len(df):
        print("Max. length of data from from start to end should be %d, but length is %d" % (len(dtimes),len(df)))
        return None
    counts = df.count(axis=0).to_numpy()
    percdata = 100*counts/len(dtimes)
    return df[ df.columns[percdata >= minperc] ]


def get_sign_level(pval):
    """Return string that corresponds to significance level based on p-value
    """
    if pval < 0.01:
        return "***"
    if pval < 0.05:
        return "** "
    if pval < 0.1:
        return "*  "
    return "   "

