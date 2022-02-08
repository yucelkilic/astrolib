# -*- coding: utf-8 -*-

import glob
import numpy as np
import pandas as pd
import paramiko
import os
import sqlite3
import folium
from geopy.point import Point
from fastkml.kml import KML
from .astronomy import FitsOps
from .astronomy import AstCalc
from datetime import datetime
from astropy.table import Table
import datetime as dt
from termcolor import colored
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.header import Header
import requests, json
import geopandas as gpd


class FileOps:

    def make_date(self, datestr, datefrmt='%Y-%m-%d'):

        return datetime.strptime(datestr.decode('ascii'), datefrmt)

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
                        sch_dict["FILTER"] = subsets.replace(
                            "'", "").strip().split(",")

                    elif "DURATION" in keyword:
                        duration = value.strip()
                        sch_dict["DURATION"] = duration.replace(
                            "'", "").strip().split(",")
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

    def print_sch_file(self, schs_path, dT=60, dF=4, dS=6, email_format=True,
                       pi_name=None, pi_surname=None, email_to=None, project_term=None):
        """
        Reads text file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @param dT: Time consumed for a single pointing
        @type float
        @param dF: Time consumed in changing the filters
        @type float
        @param dS: Time consumed in downloading the frame from CCD
        @type float
        @param email_format: Print e-mail format
        @type str
        @param pi_name: PI name
        @type float
        @param pi_surname: PI surname
        @type str
        @param email_to: Email to content
        @type string
        @param project_term: Project term
        @type string
        @return: array
        """

        sch_files = sorted(glob.glob("{0}/*.sch".format(schs_path)))
        project_proper = []
        total_exposure_time = []
        total_observation_time = []
        exposures = []
        filters = []

        filter_dict = {"U": "U",
                       "B": "B",
                       "V": "V",
                       "R": "R",
                       "I": "I",
                       "C": "C",
                       "1": "u'",
                       "2": "g'",
                       "3": "r'",
                       "4": "i'",
                       "5": "z'",
                       "H": "H-alpha",}

        sch_dict = {}
        for file_name in sch_files:
            sch_dict = self.read_sch_file(file_name)

            filter_and_durations = []
            durations = 0
            for i, filter in enumerate(sch_dict['FILTER']):
                try:
                    filter = filter_dict[filter]
                except:
                    filter = "Filtre girilmedi!"
                filter_and_durations.append(
                    "{0}({1})".format(filter, sch_dict['DURATION'][i]))
                try:
                    durations += float(sch_dict['DURATION'][i])
                except ValueError:
                    durations += 0
                try:
                    exposures.append(float(sch_dict['DURATION'][i]))
                except ValueError:
                    exposures.append(0)
                filters.append(sch_dict['FILTER'][i])
            total_exposure_time.append(durations * float(sch_dict['REPEAT']))
            total_observation_time.append(
                dT + (float(sch_dict['REPEAT']) * len(sch_dict['FILTER']) * (dF + dS)) +
                (float(sch_dict['REPEAT']) * durations))

            project_proper.append([sch_dict['SOURCE'], sch_dict['RA'], sch_dict['DEC'], " ".join(filter_and_durations),
                                   sch_dict['REPEAT']])

        pid = str(sch_dict['SOURCE']).split("_")[0]

        project_proper_np = np.asarray(project_proper)

        project_proper_tbl = Table(project_proper_np, names=('Nesne',
                                                             'RA',
                                                             'Dec',
                                                             'Filtre',
                                                             'Tekrar'))
        project_proper_tbl['Nesne'].unit = ''
        project_proper_tbl['RA'].unit = 'HH:MM:SS'
        project_proper_tbl['Dec'].unit = 'DD:MM:SS'
        project_proper_tbl['Filtre'].unit = 's'
        project_proper_tbl['Tekrar'].unit = 'cnt'

        project_proper_tbl = "".join(project_proper_tbl.pformat(html=True, max_lines=32, max_width=-1, align="<")).replace("table",
                                                                                                      'table cellspacing="0" cellpadding="6" border="1"')

        used_filters = sorted(set(filters))
        # used_filters = ''

        if '' in used_filters:
            note = "Not: Filtre/Poz süreleri belirtilmemiş nesneleriniz bulunmaktadır. " \
                   "Lütfen en kısa sürede talep ettiğiniz filtre/poz sürelerini iletiniz."
            note += "<br>"
        else:
            note = ""

        if email_format is True:
            if project_term is not None:
                if "A" in project_term:
                    project_term = project_term + " (Şubat - Mayıs)"
                elif "B" in project_term:
                    project_term = project_term + " (Haziran - Eylül)"
                elif "C" in project_term:
                    project_term = project_term + " (Ekim - Ocak)"

                message = "<strong>{pid}</strong> nolu projeniz kapsamında, <strong>{project_term}</strong> ".format(
                    pid=pid,
                    project_term=project_term)
                message += "döneminde gözlenecek olan nesneleriniz TUG PTS sisteminde aşağıdaki gibi kayıtlıdır. " \
                           "Herhangi bir düzeltme talebinizin olmaması durumunda " \
                           "gözlemleriniz bu şekilde gerçekleştirilecektir."
            else:
                message = "Talep ettiğiniz değişiklik sisteme aşağıdaki gibi uygulanmıştır."

            html = """\
                        <html>
                            <body>
                            <div style="margin-left: 1em">
                                <p>
                                    Sayın <strong>{pi_name} {pi_surname}</strong>,<br><br>
    
                                    {message}
                                    <p style="color:red">
                                        <br>
                                        Sistemde bir uyumsuzluk olduğunu düşünüyorsanız lütfen bildiriniz.
                                    </p>
                                </p>
                                Toplam nesne sayısı: <strong>{total_object_count}</strong>.
                                <br>
                                Toplam poz süresi: <strong>{total_exposure_time}</strong> saniye.
                                <br>
                                Toplam gözlem zamanı: <strong>{total_observation_time}</strong> saniye.
                                <br><br>
                                <div class="row">
                                    <div class="col-md-5">
                                        {project_table}
                                    </div>
                                </div>
                                <p style="color:red">
                                    {note}
                                </p>
                                <p>
                                    Bilgilerinize sunar,
                                    <br>
                                    İyi çalışmalar dileriz.
                                    <br><br>
                                </p>
                                TÜBİTAK Ulusal Gözlemevi (TUG) - Robotik T60 Gözlem Kontrol Sistemi (v2.0-beta) - © {year}
                                <br>
                            </div>
                    </body>
                </html>
        """.format(
                pi_name=pi_name,
                pi_surname=pi_surname,
                message=message,
                total_exposure_time=sum(total_exposure_time),
                total_observation_time=sum(total_observation_time),
                total_object_count=len(project_proper_np),
                project_table=project_proper_tbl,
                note=note,
                year=datetime.today().year)

            if email_to is not None:
                self.send_email(receivers=[email_to],
                                subject="TUG Robotik T60 Teleskobu {0} Dönemi Nesne Kontrol Talebi - {1}".format(
                                    project_term,
                                    pid),
                                content=html)

        return html

    def list_all_exposures(self, schs_path):
        """
        Reads text file into numpy array.
        @param schs_path: SCH files path
        @type schs_path: str
        @return: list
        """

        sch_files = sorted(glob.glob("{0}/*.sch".format(schs_path)))
        project_proper = []
        total_exposure_time = []
        total_observation_time = []
        exposures = []
        filters = []

        filter_dict = {"U": "U",
                       "B": "B",
                       "V": "V",
                       "R": "R",
                       "I": "I",
                       "C": "C",
                       "1": "u'",
                       "2": "g'",
                       "3": "r'",
                       "4": "i'",
                       "5": "z'",
                       "H": "H-alpha", }

        for file_name in sch_files:
            sch_dict = self.read_sch_file(file_name)

            filter_and_durations = []
            durations = 0
            for i, filter in enumerate(sch_dict['FILTER']):
                filter = filter_dict[filter]
                filter_and_durations.append(
                    "{0}({1})".format(filter, sch_dict['DURATION'][i]))
                try:
                    durations += float(sch_dict['DURATION'][i])
                except ValueError:
                    durations += 0
                try:
                    exposures.append(float(sch_dict['DURATION'][i]))
                except ValueError:
                    exposures.append(0)

                filters.append(sch_dict['FILTER'][i])

        return sorted(set(filters)), sorted(set(exposures))

    def list_pts_objects(self, project_term=None, telescope="t60",
                   key=None,
                   out_file_path="./"):
        """
        Returns all objects from TUG PTS by project term.
            Parameters
            ----------
            project_term: str
                Robotic T60 Telescope project term.
            key: str
                API key.
            out_file_path: str
                Out file path.
            Returns
            -------
            'Talon sch files'
            Example:
            -------
        """

        api_uri = "http://api.pts.tug.tubitak.gov.tr/v1/{telescope}/projects/terms/{project_term}".format(
            project_term=project_term, telescope=telescope)
        headers = {
            'Content-Type': 'application/json',
            'PTS-API-KEY': str(key)
        }
        try:
            data = requests.get(api_uri, headers=headers, timeout=1)
            if (data.status_code != 200):
                pass
        except requests.exceptions.RequestException:
            pass
        except requests.exceptions.HTTPError as e:
            return e
        except requests.exceptions.ConnectionError as e:
            return e
        except requests.exceptions.Timeout as e:
            return e

        projects = data.json()['content']

        # print(projects)

        pts_ak_score_dict = {}
        for project in projects:
            # print(project['parent_id'], float(project['ak_score']))
            if project['parent_id'] == 0:
                pid = project['id']
            else:
                pid = project['parent_id']
            print(project['objects'])

        return True

    def pts_to_sch_files(self, project_term=None,
                   key=None,
                   out_file_path="./"):
        """
        Returns all objects from TUG PTS by project term.
            Parameters
            ----------
            project_term: str
                Robotic T60 Telescope project term.
            key: str
                API key.
            out_file_path: str
                Out file path.
            Returns
            -------
            'Talon sch files'
            Example:
            -------
        """

        api_uri = "http://api.pts.tug.tubitak.gov.tr/v1/t60/projects/terms/{project_term}".format(
            project_term=project_term)
        headers = {
            'Content-Type': 'application/json',
            'PTS-API-KEY': str(key)
        }
        try:
            data = requests.get(api_uri, headers=headers, timeout=1)
            if (data.status_code != 200):
                pass
        except requests.exceptions.RequestException:
            pass
        except requests.exceptions.HTTPError as e:
            return e
        except requests.exceptions.ConnectionError as e:
            return e
        except requests.exceptions.Timeout as e:
            return e

        projects = data.json()['content']

        # print(projects)

        pts_ak_score_dict = {}
        for project in projects:
            # print(project['parent_id'], float(project['ak_score']))
            if project['parent_id'] == 0:
                pid = project['id']
            else:
                pid = project['parent_id']
            pts_ak_score_dict[pid] = float(project['ak_score'])

        pts_ak_score_dict = {k: v for k, v in sorted(pts_ak_score_dict.items(), key=lambda item: item[1], reverse=True)}

        priority_dict = {}
        counter = 0
        ak_score_check = None
        for pid, ak_score in pts_ak_score_dict.items():
            if ak_score_check != ak_score:
                counter += 1
            priority_dict[pid] = counter
            ak_score_check = ak_score

        print (pts_ak_score_dict)
        print(priority_dict)

        for project in projects:
            if project['parent_id'] == 0:
                pid = project['id']
            else:
                pid = project['parent_id']

            objects = project['objects']

            for object in objects:
                target_name = object['name'].replace(".", "_")
                target_name = target_name.replace("+", "_")
                target_name = target_name.replace("-", "_")
                target_name = target_name.replace("V*", "")
                ra = object['ra']
                dec = object['dec']

                if "-" not in dec[0]:
                    if "+" not in dec[0]:
                        dec = "+" + object['dec']

                subsets = []
                durations = []

                filter_names = object['filter_names']
                pts_durations = object['durations']
                for fw in filter_names:
                    fw_index = filter_names.index(fw)
                    if float(pts_durations[fw_index]) > 0:
                        duration = pts_durations[fw_index]

                        if "u" in fw:
                            fw = "1"
                        elif "g" in fw:
                            fw = "2"
                        elif "r" in fw:
                            fw = "3"
                        elif "i" in fw:
                            fw = "4"
                        elif "H-alpha" in fw:
                            fw = "H"
                        elif "z" in fw:
                            fw = "5"

                        subsets.append(fw)

                        durations.append(str(duration))

                priority = priority_dict[pid]
                repeat = object['repeat']

                print(pid, target_name, ra, dec, ",".join(
                    subsets), ",".join(durations), str(priority), str(repeat))

                self.create_sch_file(out_file_path=out_file_path,
                                     pid=pid,
                                     target_name=target_name,
                                     ra=ra,
                                     dec=dec,
                                     subsets=",".join(subsets),
                                     durations=",".join(durations),
                                     priority=str(priority),
                                     repeat=str(repeat))
        return True


    def csv_to_sch_files(self, csv_file, out_file_path="./", delimiter=",",
                         priority=[]):
        """
        Reads object list and creates sch_files.
        @param csv_file: Text file name and path
        @type csv_file: str
        @return: file
        """

        objects = self.read_table_as_array(csv_file, delimiter=delimiter)
        filters_in_wheel = ["C", "U", "B", "V", "R", "I",
                            "u", "g", "r", "i", "z", "H-alpha", "H-beta"]

        for object in objects:
            aco = AstCalc()

            pid = object[0]
            target_name = object[1].replace(".", "_")
            target_name = target_name.replace("+", "_")
            target_name = target_name.replace("-", "_")
            ra = aco.deg2hmsdms(object[2], object[3]).split(" ")[0][:-3]
            dec = aco.deg2hmsdms(object[2], object[3]).split(" ")[1][:-3]

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

            if len(priority) > 0:
                prior = priority.index(int(pid)) + 1
            else:
                prior = 1
            repeat = object[5]

            print(pid, target_name, ra, dec, ",".join(
                subsets), ",".join(durations), prior, repeat)

            self.create_sch_file(out_file_path=out_file_path,
                                 pid=pid,
                                 target_name=target_name,
                                 ra=ra,
                                 dec=dec,
                                 subsets=",".join(subsets),
                                 durations=",".join(durations),
                                 priority=prior,
                                 repeat=repeat)
        return True

    def get_unique_filter_duration(self, csv_file, delimiter=","):
        """
        Reads object list and returns filters and durations.
        @param csv_file: Text file name and path
        @type csv_file: str
        @return: file
        """

        objects = self.read_table_as_array(csv_file, delimiter=delimiter)
        filters_in_wheel = ["C", "U", "B", "V", "R", "I",
                            "u", "g", "r", "i", "z", "H-alpha", "H-beta"]

        subsets = []
        durations = []

        for object in objects:
            aco = AstCalc()
            for fw in filters_in_wheel:
                if fw in object[4]:
                    filters = object[4].split(";")
                    for filter_and_duration in filters:
                        if fw in filter_and_duration:
                            filter, duration = filter_and_duration.split("=")
                            filter = filter.replace("{", "")
                            duration = duration.replace("}", "")
                            if filter not in subsets:
                                subsets.append(filter)
                            if float(duration) not in durations:
                                durations.append(float(duration))

        return {'filters': subsets, 'durations': sorted(durations)}

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

            print(">>> {0}_{1}.sch file created in {2}".format(
                pid, target_name, out_file_path))

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

        filters = sch_dict['FILTER']

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

    def read_sextractor_file(self, file_name):
        """
        Reads A-Track result file into numpy array.
        @param file_name: Text file name and path
        @type file_name: str
        @return: array
        """

        try:
            with open(file_name) as fp:
                columns = []
                units = []
                for cnt, line in enumerate(fp):
                    if "#" in line:
                        line = line.replace("\n", "")
                        line = line.split()
                        columns.append(line[2])
                        if "[" in line[-1]:
                            unit = line[-1].replace("[", "")
                            unit = unit.replace("]", "")
                            units.append(unit)
                        else:
                            units.append("")
                        # print("Line {}: {}".format(cnt, line.split()))
                    else:
                        break

            data = np.genfromtxt(file_name,
                                 comments='#',
                                 invalid_raise=False,
                                 delimiter=None)

            np_result = (data[~np.isnan(data).any(axis=1)])

            tbl_result = Table(np_result,
                               names=tuple(columns))

            for k, column in enumerate(columns):
                tbl_result[column].unit = units[k]

            return tbl_result

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
                    f2db = self.fitshead_to_database(
                        fileattr.filename, sqlite_file=sqlite_file)
                    if f2db is False:
                        print(">>> Something wrong with the FITS file: {0}".format(
                            fileattr.filename))
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
                        'fwhmvs',
                        'gain',
                        'rdnoise',
                        'readtime',
                        'mjd',
                        'hjd']

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
                    band = str(fo.get_header("filter")).strip()[0]
                    try:
                        if "+" in fits_name:
                            pid = fits_name[(fits_name.index(
                                "+") + 5):(fits_name.index(band) - 3)]
                        elif "-" in fits_name:
                            pid = fits_name[(fits_name.index(
                                "-") + 5):(fits_name.index(band) - 3)]
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
                if (keyword_value is None) or (keyword_value is False):
                    keyword_value = -9999
                    
                table_headers.append(keyword)
                keyword_values.append(keyword_value)

        gain_index = table_headers.index('gain')
        hjd_index = table_headers.index('hjd')
        mjd_index = table_headers.index('mjd')
        rdnoise_index = table_headers.index('rdnoise')
        readtime_index = table_headers.index('readtime')

        gain = keyword_values[gain_index]
        hjd = keyword_values[hjd_index]
        mjd = keyword_values[mjd_index]
        rdnoise = keyword_values[rdnoise_index]
        readtime = keyword_values[readtime_index]

        table_headers.remove('gain')
        table_headers.remove('hjd')
        table_headers.remove('mjd')
        table_headers.remove('rdnoise')
        table_headers.remove('readtime')

        table_headers.append('gain')
        table_headers.append('rdnoise')
        table_headers.append('readtime')
        table_headers.append('mjd')
        table_headers.append('hjd')

        keyword_values.remove(gain)
        keyword_values.remove(hjd)
        keyword_values.remove(mjd)
        keyword_values.remove(rdnoise)
        keyword_values.remove(readtime)

        keyword_values.append(gain)
        keyword_values.append(rdnoise)
        keyword_values.append(readtime)
        keyword_values.append(mjd)
        keyword_values.append(hjd)
        
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


    def create_connection(self, db_file):
        """ create a database connection to the SQLite database
            specified by the db_file
        :param db_file: database file
        :return: Connection object or None
        """
        try:
            conn = sqlite3.connect(db_file, isolation_level=None)
            return conn
        except Error as e:
            print(e)

        return None

    def select_items_by_date(conn, database_table, start_jd, end_jd):
        """
        Query tasks by priority
        :param conn: the Connection object
        :param priority:
        :return:
        """

        cur = conn.cursor()
        cur.execute("SELECT * FROM {tbl} WHERE jd>={start_jd} AND jd<={end_jd} ORDER BY jd".format(tbl=database_table,
                                                                                                   start_jd=start_jd,
                                                                                                   end_jd=end_jd))

        rows = cur.fetchall()

        return Table(np.array(rows), names=('dateTime',
                                           'xfactor',
                                           'yfactor',
                                           'exptime',
                                           'pid',
                                           'object_name',
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
                                            'fwhmvs'))

    def srg_ephemeris_reader(self, ephemeris_file):
        """
        Reads SRG satellite  ephemeris file.
            Parameters
            ----------
            ephemeris_file: file object
                Ephemeris file.
            Returns
            -------
            'A pandas object'
            Example:
            -------
            >>> from astrolib import visuals
            >>> from astrolib import io
            >>> fo = io.FileOps()
            >>> fo.srg_ephemeris_reader("ephemeris_file.txt")
        """

        srg_ephem = np.genfromtxt(ephemeris_file,
                                  comments='=',
                                  delimiter=None,
                                  dtype="U",
                                  skip_header=5)

        srg_ephem_pd = pd.DataFrame(srg_ephem,
                                    columns=['Date',
                                             'Time',
                                             'Az',
                                             'Um',
                                             'RA2000_HH',
                                             'RA2000_MM',
                                             'RA2000_SS',
                                             'DECL2000_DD',
                                             'DECL2000_MM',
                                             'DECL2000_SS',
                                             'RARate',
                                             'DECLRate',
                                             'HourAng_HH',
                                             'HourAng_MM',
                                             'HourAng_SS',
                                             'Phase',
                                             'Illum',
                                             'SunAng',
                                             'Mag',
                                             'SDRa',
                                             'SDDecl'])

        srg_ephem_pd['Date-Time'] = pd.to_datetime(srg_ephem_pd['Date'] + ' ' + srg_ephem_pd['Time'], format='%y-%m-%d %H:%M')
        srg_ephem_pd['RA2000'] = srg_ephem_pd['RA2000_HH'] + ':' + srg_ephem_pd['RA2000_MM'] + ':' + \
                                 srg_ephem_pd['RA2000_SS']
        srg_ephem_pd['DECL2000'] = srg_ephem_pd['DECL2000_DD'] + ':' + srg_ephem_pd['DECL2000_MM'] + ':' + \
                                 srg_ephem_pd['DECL2000_SS']
        srg_ephem_pd['HourAng'] = srg_ephem_pd['HourAng_HH'] + ':' + srg_ephem_pd['HourAng_MM'] + ':' + \
                                 srg_ephem_pd['HourAng_SS']

        srg_astro_table = Table.from_pandas(srg_ephem_pd)

        return(srg_astro_table[['Date-Time',
                             'Az',
                             'Um',
                             'RA2000',
                             'DECL2000',
                             'RARate',
                             'DECLRate',
                             'HourAng',
                             'Phase',
                             'Illum',
                             'SunAng',
                             'Mag',
                             'SDRa',
                             'SDDecl']])

    def send_email(self, sender=str(Header('TUG Robotik T60 Teleskobu <tug.t60@tug.tubitak.gov.tr>')),
                   receivers=[],
                   mail_server='',
                   subject=None,
                   content=None):
        msg = MIMEMultipart('alternative')
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = ", ".join(receivers)

        html_msg = MIMEText(content, 'html')
        msg.attach(html_msg)

        try:
            smtpObj = smtplib.SMTP(mail_server, 25)
            smtpObj.sendmail(sender, receivers, msg.as_string())
            smtpObj.quit()
            print(colored("Successfully sent email to {0}".format(receivers[0]), "green"))
            return True
        except:
            print(colored("Error: unable to send email to {0}".format(receivers[0]), "red"))
            return False

    def read_kml(self, kml_file):
        """
        Reads KML files.
            Parameters
            ----------
            fname: file object
                KML file.
            Returns
            -------
            'list'
        """
        kml = KML()
        kml.from_string(open(kml_file, "rb").read())
        points = dict()
        pairs = []
        placemark_names = []

        for feature in kml.features():
            for placemark in feature.features():
                pairs = []
                for j, k, i in placemark.geometry.coords:
                    pairs.append((k, j))

                if placemark.name in placemark_names:
                    points[placemark.name + "1"] = pairs
                    placemark_names.append(placemark.name + "1")
                else:
                    points[placemark.name] = pairs
                    placemark_names.append(placemark.name)
        return points

    def read_kml_via_geopd(self, kml_file):
        """
        Reads KML files.
            Parameters
            ----------
            fname: file object
                KML file.
            Returns
            -------
            'list'
        """
        # Enable fiona driver
        gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'

        # Read file
        df = gpd.read_file(kml_file, driver='KML')

        points = dict()
        placemark_names = ['Center of shadow',
                           'Body shadow limit1',
                           'Body shadow limit2',
                           'Uncertainty1',
                           'Uncertainty2']

        for idx, feature in enumerate(df['geometry']):
            pairs = []
            for j, k, i in feature.coords:
                pairs.append([k, j])

            points[placemark_names[idx]] = pairs

        return points

    def create_occultation_map(self,
                               location_file,
                               kml_file,
                               sep="|", header=0,
                               observatory_column_keyword="Observatory",
                               latitude_column_keyword="Latitude",
                               longitude_column_keyword="Longitude", save_map=False):
        """
        Creates occultation map.
            Parameters
            ----------
            location_file: file object
                Location file.
            sep: str
                Seperator of row data.
            header: int
                Row position of header for input file.
            latitude_column_keyword: str
                Keyword of latitude column.
            longitude_column_keyword: str
                Keyword of longitude column.
            Returns
            -------
            'A map object'
            Example:
            -------
            >>> from astrolib import io
            >>> fo = io.FileOps()
            >>> fo.create_occultation_map(kml_file="2002KX14_20200526_NIMAv7_LuckyStar.kmz", location_file="2002KX14_locations.txt")
        """
        for_map = pd.read_csv(location_file, sep=sep, header=header)
        locations = self.read_kml(kml_file)

        middle_point = int(len(locations['Center of shadow']) / 2)

        # Make an empty map
        m = folium.Map(location=(40, 35), zoom_start=5)


        # I can add marker one by one on the map
        for i in range(0, len(for_map)):
            long = str(for_map[longitude_column_keyword][i]).strip()
            lat = str(for_map[latitude_column_keyword][i]).strip()
            observatory = str(for_map[observatory_column_keyword][i]).strip()

            long = long.replace("°", " ")

            lat = lat.replace("°", " ")

            p = Point('''{lat} {long}'''.format(lat=lat, long=long))
            # print(observatory, p.format_decimal)

            latitude, longitude, altitude = p

            folium.Marker([latitude, longitude], popup=observatory).add_to(m)

        folium.PolyLine(locations['Body shadow limit'], popup="Body shadow upper limit").add_to(m)
        folium.PolyLine(locations['Body shadow limit1'], popup="Body shadow bottom limit").add_to(m)
        folium.PolyLine(locations['Center of shadow'][1:], popup="Center of shadow", color='green').add_to(m)
        folium.PolyLine(locations['Uncertainty'], popup="Uncertainty", color='red', dash_array='10').add_to(m)
        folium.PolyLine(locations['Uncertainty1'], popup="Uncertainty ", color='red', dash_array='10').add_to(m)

        if save_map:
            path_name, ext = os.path.splitext(location_file)
            m.save('{}.html'.format(path_name))

        return m

    def occultation_path(self, kml_file):
        """
        Creates occultation map.
            Parameters
            ----------
            location_file: file object
                Location file.
            sep: str
                Seperator of row data.
            header: int
                Row position of header for input file.
            latitude_column_keyword: str
                Keyword of latitude column.
            longitude_column_keyword: str
                Keyword of longitude column.
            Returns
            -------
            'A map object'
            Example:
            -------
            >>> from astrolib import io
            >>> fo = io.FileOps()
            >>> fo.create_occultation_map(kml_file="2002KX14_20200526_NIMAv7_LuckyStar.kmz")
        """
        locations = self.read_kml_via_geopd(kml_file)

        return {'body_upper_limit': locations['Body shadow limit1'],
                'body_bottom_limit': locations['Body shadow limit2'],
                'body_center': locations['Center of shadow'],
                'upper_uncertain':  locations['Uncertainty1'],
                'bottom_uncertain':  locations['Uncertainty2']}
