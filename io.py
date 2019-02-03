# -*- coding: utf-8 -*-

import glob
import numpy as np
import paramiko
import os
import sqlite3
from .astronomy import FitsOps
from .astronomy import AstCalc
from datetime import datetime


class FileOps:

    def make_date(self, datestr, datefrmt='%Y-%m-%d'):

        return(datetime.strptime(datestr.decode('ascii'), datefrmt))

    def read_file_as_array(self, file_name):

        """
        Reads text file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:
            return (np.genfromtxt(file_name,
                                  comments='#',
                                  delimiter=' | ',
                                  dtype="U"))
        except Exception as e:
            print(e)

    def read_table_as_array(self, file_name, delimiter=None):
        """
        Reads A-Track result file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        data = np.genfromtxt(file_name,
                             delimiter=delimiter,
                             dtype="|U30")
        return data



    def read_sch_file(self, file_name):

        """
        Reads text file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:
            with open(file_name) as fin:
                sch_dict = {}
                for line in fin:
                    keyword_value = line.replace("\n", "")
                    if "=" in keyword_value:
                        keyword, value = keyword_value.split("=")
                    else:
                        continue

                    if "TITLE" in keyword:
                        sch_dict["TITLE"] = value.strip()
                    elif "OBSERVER" in keyword:
                        sch_dict["OBSERVER"] = value.strip()
                    elif "BLOCK" in keyword:
                        sch_dict["BLOCK"] = value.strip()
                    elif "SOURCE" in keyword:
                        sch_dict["SOURCE"] = value.strip()
                    elif " RA" in keyword:
                        sch_dict["RA"] = value.strip()
                    elif "DEC" in keyword:
                        sch_dict["DEC"] = value.strip()
                    elif "EPOCH" in keyword:
                        sch_dict["EPOCH"] = value.strip()
                    elif "FILTER" in keyword:
                        subsets = value.strip()
                        sch_dict["FILTER"] = subsets.replace("'", "").strip().split(",")
                        
                    elif "DURATION" in keyword:
                        duration = value.strip()
                        sch_dict["DURATION"] = duration.replace("'", "").strip().split(",")
                    elif "BINNING" in keyword:
                        sch_dict["BINNING"] = value.strip()
                    elif "SUBIMAGE" in keyword:
                        sch_dict["SUBIMAGE"] = value.strip()
                    elif "LSTDELTA" in keyword:
                        sch_dict["LSTDELTA"] = value.strip()
                    elif "PRIORITY" in keyword:
                        sch_dict["PRIORITY"] = value.strip()
                    elif "COMPRESS" in keyword:
                        sch_dict["COMPRESS"] = value.strip()
                    elif "IMAGEDIR" in keyword:
                        sch_dict["IMAGEDIR"] = value.strip()
                    elif " REPEAT" in keyword:
                        sch_dict["REPEAT"] = value.strip()
                    elif "BLOCKREPEAT" in keyword:
                        sch_dict["BLOCKREPEAT"] = value.strip()
                        
                return sch_dict
        except Exception as e:
            print(e)
            return False

    def csv_to_sch_files(self, csv_file, out_file_path="./", delimiter=","):

        """
        Reads object list and creates sch_files.
        @param csv_file: Text file name and path
        @type csv_file: str
        @return: file
        """

        objects = self.read_table_as_array(csv_file, delimiter=",")
        filters_in_wheel = ["C", "U", "B", "V", "R", "I",
                            "u", "g", "r", "z", "H-alpha", "H-beta"]

        for object in objects:
            aco = AstCalc()

            pid = object[0]
            target_name = object[1].replace(".", "_")
            ra = aco.deg2hmsdms(object[2], object[3]).split(" ")[0]
            dec = aco.deg2hmsdms(object[2], object[3]).split(" ")[1]

            subsets = []
            durations = []
            for fw in filters_in_wheel:
                if fw in object[4]:
                    filters = object[4].split(";")
                    for filter_and_duration in filters:
                        if fw in filter_and_duration:
                            filter, duration = filter_and_duration.split("=")
                            filter = filter.replace("{", "")
                            duration = duration.replace("}", "")
                            subsets.append(filter)
                            durations.append(duration)

            priority = 1
            repeat = object[5]

            print(pid, target_name, ra, dec, ",".join(subsets), ",".join(durations), priority, repeat, sep=",")

            self.create_sch_file(out_file_path=out_file_path,
                                 pid=pid,
                                 target_name=target_name,
                                 ra=ra,
                                 dec=dec,
                                 subsets=",".join(subsets),
                                 durations=",".join(durations),
                                 priority=1,
                                 repeat=repeat)
        return True


    def create_sch_file(self, out_file_path="./", pid="3141", target_name="tug", ra=None,
                        dec=None,
                        subsets="U,B,V,R,I",
                        durations="60,60,60,60,60",
                        priority=1,
                        repeat=1):

        """
        Reads A-Track result file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:

            sch_template = """! Generated by TUG/astrolib

TITLE = {0}_{1}
OBSERVER = 'TUG'
        
BLOCK = '1'
    SOURCE = {0}_{1}
    RA = {2}
    DEC = {3}
    EPOCH = 2000
    FILTER = '{4}'
    DURATION = '{5}'
    BINNING = '1,1'
    SUBIMAGE = '0,0,2048,2048'
    LSTDELTA = '3'
    PRIORITY = {6}
    COMPRESS = 0
    IMAGEDIR = '/usr/local/telescope/user/images'
    REPEAT = {7}
    /
BLOCKREPEAT = 1
""".format(pid,
                       target_name,
                       ra,
                       dec,
                       subsets,
                       durations,
                       priority,
                       repeat)


            with open("{0}/{1}_{2}.sch".format(out_file_path,
                                               pid, target_name), "w") as file:
                file.write(sch_template)
                file.close()

            print(">>> {0}_{1}.sch file created in {2}".format(pid, target_name, out_file_path))

            return True
            
        except Exception as e:
            print(e)
            return False

    def sch_to_database(self, sch_file,
                             sqlite_file="db.sqlite3",
                             table_name="portal_schedule",
                             keywords=None):

        print(">>> SCH2DB: {0}".format(sch_file))

        if keywords is None:
            keywords = ['source_name',
                        'ra',
                        'dec',
                        'epoch',
                        'filterU',
                        'filterB',
                        'filterV',
                        'filterR',
                        'filterI',
                        'filterC',
                        'priority',
                        'repeat',
                        'is_ongoing_project',
                        'project_id']

        sch_dict = self.read_sch_file(sch_file)

        filters =sch_dict['FILTER']

        try:
            u_filter = sch_dict['DURATION'][filters.index("U")]
        except ValueError:
            u_filter = 0

        try:
            b_filter = sch_dict['DURATION'][filters.index("B")]
        except ValueError:
            b_filter = 0

        try:
            v_filter = sch_dict['DURATION'][filters.index("V")]
        except ValueError:
            v_filter = 0

        try:
            r_filter = sch_dict['DURATION'][filters.index("R")]
        except ValueError:
            r_filter = 0

        try:
            i_filter = sch_dict['DURATION'][filters.index("I")]
        except ValueError:
            i_filter = 0

        try:
            c_filter = sch_dict['DURATION'][filters.index("C")]
        except ValueError:
            c_filter = 0

        pid = sch_dict['SOURCE'].split("_")[0]

        print(["_".join(sch_dict['SOURCE'].split("_")[1:]),
               sch_dict['RA'],
               sch_dict['DEC'],
               sch_dict['EPOCH'],
               u_filter,
               b_filter,
               v_filter,
               r_filter,
               i_filter,
               c_filter,
               sch_dict['PRIORITY'],
               sch_dict['REPEAT'],
               1,
               pid])

        # Connecting to the database file
        conn = sqlite3.connect(sqlite_file)
        c = conn.cursor()

        c.execute("INSERT OR IGNORE INTO {tn} {cn} VALUES {vls}".format(
            tn=table_name,
            cn=tuple(keywords),
            vls=tuple(["_".join(sch_dict['SOURCE'].split("_")[1:]),
                       sch_dict['RA'],
                       sch_dict['DEC'],
                       sch_dict['EPOCH'],
                       u_filter,
                       b_filter,
                       v_filter,
                       r_filter,
                       i_filter,
                       c_filter,
                       sch_dict['PRIORITY'],
                       sch_dict['REPEAT'],
                       1,
                       pid])))

        conn.commit()
        conn.close()
        return (True)

    def read_res(self, file_name):

        """
        Reads A-Track result file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:
            data = np.genfromtxt(file_name,
                                 comments='#',
                                 skip_header=2,
                                 invalid_raise=False,
                                 delimiter=None,
                                 usecols=(0, 1, 3, 4, 5))
            return (data[~np.isnan(data).any(axis=1)])
        except Exception as e:
            print(e)

    def get_file_list(self, dir_name):

        """
        List FITS images in a folder into a numpy array.
        @param dir_name: Directory of FITS images
        @type dir_name: str
        @return: array
        """

        images = sorted(glob.glob(dir_name + '/*.fit*'))
        return (images)

    def get_fits_from_server(self,
                             hostname,
                             username,
                             password,
                             dirname="/mnt/data/images",
                             fits_ext=".fts",
                             header2sqlite=False,
                             sqlite_file="gozlemler.sqlite"):
        
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        try:
            ssh.connect(hostname=hostname,
                        username=username,
                        password=password)
        except paramiko.SSHException:
            print("Connection Failed")
            quit()

        sftp = ssh.open_sftp()
        try:
            sftp.chdir(dirname)
        except FileNotFoundError:
            print("Folder not found!")

        ret = False

        for fileattr in sftp.listdir_attr():
            if not os.path.exists(fileattr.filename) and \
                    fits_ext in fileattr.filename:
                sftp.get(fileattr.filename, fileattr.filename)
                if header2sqlite is True:
                    f2db = self.fitshead_to_database(fileattr.filename, sqlite_file=sqlite_file)
                    if f2db is False:
                        print(">>> Something wrong with the FITS file: {0}".format(fileattr.filename))
                        continue
                ret = True
                print("{0} => {1}".format(fileattr.filename,
                                          fileattr.filename))

        print("Done")
        ssh.close()
        if ret is False:
            print("No file(s) found!")

        return(ret)

    def get_latest_fits_from_server(self,
                                    hostname,
                                    username,
                                    password,
                                    dirname="/mnt/data/images",
                                    fits_ext=".fts",
                                    header2sqlite=False,
                                    sqlite_file="observation.db"):

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        try:
            ssh.connect(hostname=hostname,
                        username=username,
                        password=password)
        except paramiko.SSHException:
            print("Connection Failed")
            quit()

        sftp = ssh.open_sftp()
        try:
            sftp.chdir(dirname)
        except FileNotFoundError:
            print("No such folder!")
            ssh.close()
            raise SystemExit


        ret = False
        latest = 0
        latestfile = None

        for fileattr in sftp.listdir_attr():
            if fits_ext in fileattr.filename and fileattr.st_mtime > latest:
                latest = fileattr.st_mtime
                latestfile = fileattr.filename

        if latestfile is not None and \
                os.path.exists(fileattr.filename) is False:
            sftp.get(latestfile, latestfile)
            print("{0} => {1}".format(latestfile,
                                      latestfile))
            ret = True

        if header2sqlite is True:
            self.fitshead_to_database(latestfile, sqlite_file=sqlite_file)

        print("Done")
        ssh.close()
        if ret is False:
            print("No file(s) found!")

        return (latestfile)

    def fitshead_to_database(self, fits_file,
                             sqlite_file="observations.db",
                             table_name="T60",
                             keywords=None):

        print(">>> FITS2DB: {0}".format(fits_file))

        if keywords is None:
            keywords = ['xfactor',
                        'yfactor',
                        'exptime',
                        'object',
                        'priority',
                        'instrume',
                        'jd',
                        'date-obs',
                        'time-obs',
                        'lst',
                        'latitude',
                        'elevatio',
                        'azimuth',
                        'ha',
                        'ra',
                        'dec',
                        'objra',
                        'objdec',
                        'epoch',
                        'equinox',
                        'filter',
                        'camtemp',
                        'focuspos',
                        'wxtemp',
                        'wxpres',
                        'wxwndspd',
                        'wxwnddir',
                        'wxhumid',
                        'biascor',
                        'thermcor',
                        'flatcor',
                        'badpxcor',
                        'fwhmh',
                        'fwhmhs',
                        'fwhmv',
                        'fwhmvs']

        # Connecting to the database file
        conn = sqlite3.connect(sqlite_file)
        c = conn.cursor()

        try:
            fo = FitsOps(fits_file)
        except Exception as e:
            print(e)
            return False
            
        keyword_values = []
        table_headers = []
        fits_name = os.path.basename(fits_file)
        keyword_values.append(fits_name)
        table_headers.append("fits_name")
        for keyword in keywords:
            if keyword == "object":
                value = fo.get_header(keyword)
                if value is not None:
                    object_name = value
                    band = str(fo.get_header("filter")).strip()
                    try:
                        if "+" in fits_name:
                            pid = fits_name[(fits_name.index("+") + 5):(fits_name.index(band) - 3)]
                        elif "-" in fits_name:
                            pid = fits_name[(fits_name.index("-") + 5):(fits_name.index(band) - 3)]
                        elif "ldt" in fits_name:
                            pid = "STD"
                        else:
                            pid = -9999
                    except ValueError:
                        pid = -9999

                    print("PID: {0}".format(pid))
                else:
                    pid = -9999
                    object_name = -9999
                table_headers.append("pid")
                table_headers.append("object_name")
                keyword_values.append(pid)
                keyword_values.append(object_name)
            else:
                keyword_value = fo.get_header(keyword)
                if keyword_value is None:
                    keyword_value = -9999
                table_headers.append(keyword)
                keyword_values.append(keyword_value)

        print(table_headers)
        print(keyword_values)
                            
        c.execute("INSERT OR IGNORE INTO {tn} {cn} VALUES {vls}".format(
            tn=table_name,
            cn=tuple(table_headers),
            vls=tuple(keyword_values)))

        conn.commit()
        conn.close()

        return(True)

    def find_if_in_database_id(self, database, idd):
        """
        Search detected asteroids ID in the MPCORB.DAT database for MPC report.
        @param database: MPCORB.DAT path
        @type database: str
        @param idd: Asteroid's ID
        @type idd: str
        @return: str
        """

        ret = ""
        try:
            f = open(database, "r")
            for i in f:
                ln = i.replace("\n", "").split()
                try:
                    if "({0})".format(idd) == ln[21]:
                        id_name = ln[0]
                        if len(id_name) > 5:
                            ret = "     " + id_name
                        else:
                            ret = id_name
                except:
                    continue
            f.close()
        except Exception as e:
            print(e)

        return (ret)

    def find_if_in_database_name(self, database, name):

        """
        Search detected asteroids ID by the name in the
        MPCORB.DAT database for MPC report.

        @param database: MPCORB.DAT path
        @type database: str
        @param name: Asteroid's name
        @type name: str
        @return: str
        """

        ret = ""
        try:
            f = open(database, "r")
            for i in f:
                ln = i.replace("\n", "").split()
                try:
                    if len(ln[23]) < 8:
                        combname = "{0} {1}".format(ln[22], ln[23])
                    else:
                        combname = ln[22]

                    if name == combname:
                        id_name = ln[0]
                        if len(id_name) > 5:
                            ret = "     " + id_name
                        else:
                            ret = id_name
                except:
                    continue
            f.close()
        except Exception as e:
            print(e)

        return (ret)
