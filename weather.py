from datetime import timedelta
from datetime import date
import time

import numpy as np
import urllib
import re

from astropy.table import Table, vstack
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
from astroplan import Observer
import sqlite3
from sqlite3 import Error


class Weather:

    def daterange(self, start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(n)

    def get_current_meteo_data(self, station="RTT150"):

        if station.upper() == "T60":
            meteo_feed = feedparser.parse("http://t60meteo.tepe.tug.tubitak.gov.tr/wxrss.xml")
        elif station.upper() == "RTT150":
            meteo_feed = feedparser.parse("http://rtt150meteo.tepe.tug.tubitak.gov.tr/wxrss.xml")
        elif station.upper() == "T100":
            meteo_feed = feedparser.parse("http://t100meteo.tug.tubitak.gov.tr/wxrss.xml")
        else:
            print("No station has found!")
            raise SystemExit

        for entry in meteo_feed.entries:

            article_title = entry.title
            article_published_at = entry.published  # Unicode string

            content = entry.content

            print("{}".format(article_title))
            print("Published at {}".format(article_published_at))

            for line in content[0].value.split("<br />"):
                if ("<p>" or "</p>" or "\n") not in line:
                    line = line.replace("</p>", "")

                    if "Wind Chill:" in line:
                        windchil = line.split(" ")[-2]

                    if "Heat Index:" in line:
                        heatindex = line.split(" ")[-2]

                    if "Humidity:" in line:
                        humidity = line.split(" ")[-2]

                    if "Dewpoint:" in line:
                        dewpoint = line.split(" ")[-2]

                    if "Barometer:" in line:
                        barometer = line.split(" ")[-2]

                    if "Wind:" in line:
                        winddir = line.split(" ")[-4]
                        windspeed = line.split(" ")[-2]

                    if "Rain Today:" in line:
                        raintoday = line.split(" ")[-2]

                    if "Rain Rate:" in line:
                        rainrate = line.split(" ")[-2]

        print("Windspeed: (km/h) {} at {}".format(windspeed, winddir))
        print("Humidity (%):", humidity)

        return {"Time": article_published_at,
                "Windchil": float(windchil),
                "Heatindex": float(heatindex),
                "Humidity": float(humidity),
                "Dewpoint": float(dewpoint),
                "Barometer": float(barometer),
                "Winddir": winddir,
                "Windspeed": float(windspeed),
                "Rain Today": float(raintoday),
                "Rain Rate": float(rainrate)
                }

    def read_davis_data_from_archive(self, date_obs, station="T60", db_file=None):
        """
        It reads wxview result data from meteo station.
        @param date_obs
        @type date_obs: date
        @param station: Station name: T60|T100|RTT150.
        @type station: string
        @param station: wview sqlite db file.
        @type station: file
        @return: astropy.table
        """

        year, month, day = re.split('[- :/.]', date_obs)

        if station.upper() == "T60":
            urlink = "http://t60meteo.tepe.tug.tubitak.gov.tr/index.html/" \
                "Archive/ARC-{0}-{1}-{2}.txt".format(year, month, day)
        elif station.upper() == "RTT150":
            urlink = "http://rtt150meteo.tepe.tug.tubitak.gov.tr/" \
                "ARC-{0}-{1}-{2}.txt".format(year, month, day)
        elif station.upper() == "T100":
            urlink = "http://t100meteo.tepe.tug.tubitak.gov.tr/index.html/" \
                "Archive/ARC-{0}-{1}-{2}.txt".format(year, month, day)
        else:
            print("No station has found!")
            raise SystemExit

        # print("Retriving data from: {0}".format(urlink))
        try:
            if db_file is None:
                raw_data = urllib.request.urlopen(urlink)
                dataset = np.genfromtxt(raw_data,
                                        delimiter=None,
                                        comments="---",
                                        skip_header=1,
                                        dtype="|U30")
        except:
            return False

        if db_file is None:
            # date and time merging
            date_dataset = dataset.astype(object)
            date_dataset[:, 0] += "T" + date_dataset[:, 1]

            """
            convert first column time format because davis data
            not formatted properly
            """

            ts = [time.strftime('%Y-%m-%dT%H:%M',
                                time.strptime(s[0].replace("T24:00", "T23:59"),
                                              '%Y%m%dT%H:%M')) for s in
                  date_dataset[:, :1]]

            # convert np arrar for merging rest of data
            formatted_date_column = np.asarray(ts).reshape(len(ts), 1)

            np_davis_data = np.hstack((formatted_date_column,
                                       dataset[:, 2:]))

            tbl_davis_data = Table(np_davis_data,
                                   names=('Timestamp',
                                          'Temp',
                                          'Chill',
                                          'HIindex',
                                          'Humid',
                                          'Dewpt',
                                          'Wind',
                                          'HiWind',
                                          'WindDir',
                                          'Rain',
                                          'Barom',
                                          'Solar',
                                          'ET',
                                          'UV'),
                                   dtype=('datetime64',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8',
                                          'f8'))
        else:
            conn = sqlite3.connect(sqlite_file)
            c = conn.cursor()


        return(tbl_davis_data)

    def astronomical_twilight(self, date_obs,
                              site_longitude=30.335555,
                              site_latitude=36.824166,
                              site_elevation=2500,
                              site_name="tug",
                              time_zone="Europe/Istanbul",
                              which="next"):
        """
        It calculates astronomical twilight.
        @param date_obs
        @type date_obs: date
        @param site_longitude: Site longitude.
        @type site_longitude: float
        @param site_latitude: Site latitude.
        @type site_latitude: float
        @param site_elevation: Site elevation.
        @type site_elevation: float
        @param site_name: Station name: T60|T100|RTT150.
        @type site_name: string
        @param time_zone: Time zone.
        @type time_zone: int
        @param which: Which.
        @type which: string
        @return: tuple
        """

        # TUG's location info settings
        tug = Observer(longitude=site_longitude*u.deg,
                       latitude=site_latitude*u.deg,
                       elevation=site_elevation*u.m,
                       name=site_name,
                       timezone=time_zone)

        # convert date astropy date format
        astropy_time = Time(date_obs)

        # evening tw calculate
        et = tug.twilight_evening_astronomical(astropy_time,
                                               which=which)

        # morning tw calculate
        mt = tug.twilight_morning_astronomical(astropy_time,
                                               which=which)

        # localized time conversion
        return(tug.astropy_time_to_datetime(mt),
               tug.astropy_time_to_datetime(et))

    def nautical_twilight(self, date_obs,
                              site_longitude=30.335555,
                              site_latitude=36.824166,
                              site_elevation=2500,
                              site_name="tug",
                              time_zone="Europe/Istanbul",
                              which="next"):

        # TUG's location info settings
        tug = Observer(longitude=site_longitude*u.deg,
                       latitude=site_latitude*u.deg,
                       elevation=site_elevation*u.m,
                       name=site_name,
                       timezone=time_zone)

        # convert date astropy date format
        astropy_time = Time(date_obs)

        # evening tw calculate
        et = tug.twilight_evening_nautical(astropy_time,
                                               which=which)

        # morning tw calculate
        mt = tug.twilight_morning_nautical(astropy_time,
                                               which=which)

        # localized time conversion
        return(tug.astropy_time_to_datetime(mt),
               tug.astropy_time_to_datetime(et))

    def daily_bad_weather_report(self, date_obs, station="T60",
                                 site_longitude=30.335555,
                                 site_latitude=36.824166,
                                 site_elevation=2500,
                                 site_name="tug",
                                 time_zone="Europe/Istanbul",
                                 humidity_limit=85,
                                 wind_limit=36,
                                 twilight="nautical"):

        year, month, day = re.split('[- :/.]', date_obs)

        date_obs = "{0}-{1}-{2}".format(year, month, day)

        # d = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S")
        # date_obs = d.strftime("%Y-%m-%d")


        if twilight is "nautical":
            # date from after mid

            mt_before, et_before = self.nautical_twilight(
                date_obs,
                site_longitude,
                site_latitude,
                site_elevation,
                site_name,
                time_zone)

            end_date = date(int(year), int(month), int(day)) + timedelta(1)
            mt_after, et_after = self.nautical_twilight(
                end_date.strftime("%Y-%m-%d"),
                site_longitude,
                site_latitude,
                site_elevation,
                site_name,
                time_zone)
        elif twilight is "astronomical":
            # date from after mid

            mt_before, et_before = self.astronomical_twilight(
                date_obs,
                site_longitude,
                site_latitude,
                site_elevation,
                site_name,
                time_zone)

            end_date = date(int(year), int(month), int(day)) + timedelta(1)
            mt_after, et_after = self.astronomical_twilight(
                end_date.strftime("%Y-%m-%d"),
                site_longitude,
                site_latitude,
                site_elevation,
                site_name,
                time_zone)

        davis_data_before_mid = self.read_davis_data_from_archive(
            date_obs,
            station=station)

        davis_data_after_mid = self.read_davis_data_from_archive(
            end_date.strftime("%Y-%m-%d"),
            station=station)

        if (davis_data_before_mid and davis_data_after_mid) is not False:
            davis_data = vstack([davis_data_before_mid, davis_data_after_mid])
        else:
            davis_data = False

        # calculate dark night time
        tw_dates = [et_before, mt_after]
        t = Time(tw_dates) + TimeDelta(time.timezone * -1, format='sec')
        dark_hours = (t[1].jd - t[0].jd) * 24

        if davis_data is not False:
            davis_data['jd'] = Time(davis_data['Timestamp'].astype("str")).jd

            date_obs_davis_data_jd = davis_data[(davis_data['jd'] >= t[0].jd) & (davis_data['jd'] <= t[1].jd)]
            bad_data = date_obs_davis_data_jd[(date_obs_davis_data_jd['Humid'] > humidity_limit) |
                                              (date_obs_davis_data_jd['Wind'] >= wind_limit)]

            bad_data_diff_sec = np.diff(bad_data['jd']) * 86400.0
            daily_bad_weather_times = bad_data_diff_sec[bad_data_diff_sec <= 1800.0]
            interrupt_cnt = len(bad_data_diff_sec[bad_data_diff_sec > 1800.0])

            if len(bad_data['jd']) > 0:
                first_data = 5
            else:
                first_data = 0

            total_bad_weather_time = (np.sum(daily_bad_weather_times) / 3600.0) + (interrupt_cnt * 5.0 / 60.0) \
                                     + (first_data / 60.0)

            humid_data = bad_data[(bad_data['Humid'] > humidity_limit)]
            wind_data = bad_data[(bad_data['Wind'] >= wind_limit)]

            humid_cnt = len(humid_data)
            wind_cnt = len(wind_data)

            if (humid_cnt + wind_cnt) == 1:
                weather_alert_time = 0.5
                total_bad_weather_time = (5.0 * (humid_cnt + wind_cnt) / 60.0) + weather_alert_time

            if wind_cnt == 0 and humid_cnt == 0:
                weather_status = "Clear"
            elif humid_cnt > 0 and wind_cnt == 0:
                weather_status = "Humid"
            elif wind_cnt > 0 and humid_cnt == 0:
                weather_status = "Windy"
            elif wind_cnt > 0 and humid_cnt > 0:
                weather_status = "HumidWindy"
        else:
            total_bad_weather_time = dark_hours
            weather_status = "NoConnection"

        return (date_obs,
                '{0:.2f}'.format(dark_hours),
                '{0:.2f}'.format(total_bad_weather_time),
                weather_status)

    def long_term_bad_weather_report(self, start_date,
                                     end_date,
                                     station="T60",
                                     site_longitude=30.335555,
                                     site_latitude=36.824166,
                                     site_elevation=2500,
                                     site_name="tug",
                                     time_zone="Europe/Istanbul",
                                     humidity_limit=85,
                                     wind_limit=36):
        syear, smonth, sday = re.split('[- :/.]', start_date)

        start_date = "{0}-{1}-{2}".format(syear, smonth, sday)

        eyear, emonth, eday = re.split('[- :/.]', end_date)

        end_date = "{0}-{1}-{2}".format(eyear, emonth, eday)

        start_date = date(int(syear), int(smonth), int(sday))
        end_date = date(int(eyear), int(emonth), int(eday))

        print("{0} {1} {2}\n".format("date", "dark_hours",
                                     "bad_weather_hours"))

        with open("{0}_bad_weather_report-{1}_{2}.txt".format(
                station,
                start_date,
                end_date), "a") as out:

            out.write("{0} {1} {2}\n".format("date", "dark_hours",
                                             "bad_weather_hours"))

            for single_date in self.daterange(start_date, end_date):
                odate = single_date.strftime("%Y-%m-%d")
                dt, dark_hours, bad_time, weather_status = self.daily_bad_weather_report(
                    odate,
                    station,
                    site_longitude,
                    site_latitude,
                    site_elevation,
                    site_name,
                    time_zone,
                    humidity_limit,
                    wind_limit)

                print(dt, dark_hours, bad_time, weather_status)
                out.write("{0} {1} {2}\n".format(dt,
                                                 dark_hours,
                                                 bad_time, weather_status))

        print("Long term bad weather report has been written "
              "to {0}_bad_weather_report-{1}_{2}.txt".format(
                  station,
                  start_date,
                  end_date))

        return(True)
