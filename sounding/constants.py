"""
Shared constants: station lists, model dictionaries, data source catalogue.
"""

STATIONS = {
    "OUN": ("Norman, OK", 35.22, -97.46),
    "FWD": ("Fort Worth, TX", 32.83, -97.30),
    "AMA": ("Amarillo, TX", 35.23, -101.71),
    "DDC": ("Dodge City, KS", 37.77, -99.97),
    "TOP": ("Topeka, KS", 39.07, -95.62),
    "SGF": ("Springfield, MO", 37.23, -93.40),
    "LZK": ("Little Rock, AR", 34.83, -92.26),
    "SHV": ("Shreveport, LA", 32.45, -93.84),
    "BMX": ("Birmingham, AL", 33.18, -86.77),
    "JAN": ("Jackson, MS", 32.32, -90.08),
    "LCH": ("Lake Charles, LA", 30.12, -93.23),
    "BNA": ("Nashville, TN", 36.25, -86.56),
    "ILN": ("Wilmington, OH", 39.42, -83.82),
    "DTX": ("Detroit, MI", 42.70, -83.47),
    "DVN": ("Davenport, IA", 41.61, -90.58),
    "OAX": ("Omaha, NE", 41.32, -96.37),
    "LBF": ("North Platte, NE", 41.13, -100.68),
    "ABR": ("Aberdeen, SD", 45.45, -98.41),
    "UNR": ("Rapid City, SD", 44.07, -103.21),
    "BIS": ("Bismarck, ND", 46.77, -100.75),
    "GJT": ("Grand Junction, CO", 39.12, -108.53),
    "DNR": ("Denver, CO", 39.77, -104.88),
    "ABQ": ("Albuquerque, NM", 35.04, -106.62),
    "EPZ": ("El Paso, TX", 31.87, -106.70),
    "TFX": ("Great Falls, MT", 47.46, -111.38),
    "SLC": ("Salt Lake City, UT", 40.77, -111.97),
    "BOI": ("Boise, ID", 43.57, -116.22),
    "MFR": ("Medford, OR", 42.37, -122.87),
    "OTX": ("Spokane, WA", 47.68, -117.63),
    "UIL": ("Quillayute, WA", 47.95, -124.55),
    "REV": ("Reno, NV", 39.57, -119.80),
    "VBG": ("Vandenberg, CA", 34.75, -120.57),
    "NKX": ("San Diego, CA", 32.87, -117.15),
    "TUS": ("Tucson, AZ", 32.23, -110.95),
    "FGZ": ("Flagstaff, AZ", 35.23, -111.82),
    "IAD": ("Sterling, VA", 38.98, -77.48),
    "WAL": ("Wallops Island, VA", 37.94, -75.47),
    "MHX": ("Morehead City, NC", 34.78, -76.88),
    "GSO": ("Greensboro, NC", 36.10, -79.95),
    "CHS": ("Charleston, SC", 32.90, -80.03),
    "JAX": ("Jacksonville, FL", 30.50, -81.70),
    "TBW": ("Tampa Bay, FL", 27.70, -82.40),
    "MFL": ("Miami, FL", 25.75, -80.38),
    "TLH": ("Tallahassee, FL", 30.40, -84.35),
    "BUF": ("Buffalo, NY", 42.93, -78.73),
    "ALB": ("Albany, NY", 42.75, -73.80),
    "OKX": ("Upton, NY", 40.87, -72.87),
    "GYX": ("Gray, ME", 43.90, -70.25),
    "CHH": ("Chatham, MA", 41.67, -69.97),
    "CAR": ("Caribou, ME", 46.87, -68.02),
    "PIT": ("Pittsburgh, PA", 40.53, -80.23),
    "RNK": ("Blacksburg, VA", 37.20, -80.40),
    "MPX": ("Minneapolis, MN", 44.85, -93.57),
    "GRB": ("Green Bay, WI", 44.48, -88.13),
    "ILX": ("Lincoln, IL", 40.15, -89.34),
    "APX": ("Gaylord, MI", 44.90, -84.72),
    "INL": ("International Falls, MN", 48.57, -93.38),
    "CRP": ("Corpus Christi, TX", 27.77, -97.50),
    "MAF": ("Midland, TX", 31.95, -102.18),
    "DRT": ("Del Rio, TX", 29.37, -100.92),
    "BRO": ("Brownsville, TX", 25.92, -97.42),
    "KEY": ("Key West, FL", 24.55, -81.75),
    "RIW": ("Riverton, WY", 43.07, -108.48),
    "LMN": ("Lamont, OK", 36.60, -97.48),
    "FFC": ("Peachtree City, GA", 33.36, -84.57),
    "LIX": ("Slidell, LA", 30.34, -89.83),
    "VEF": ("Las Vegas, NV", 36.05, -115.18),
    "OAK": ("Oakland, CA", 37.75, -122.22),
    "GGW": ("Glasgow, MT", 48.21, -106.63),
    "ANC": ("Anchorage, AK", 61.17, -150.02),
    "FAI": ("Fairbanks, AK", 64.82, -147.87),
    "LIH": ("Lihue, HI", 21.98, -159.35),
    "ITO": ("Hilo, HI", 19.72, -155.07),
}

# University of Wyoming station IDs (3-letter -> WMO number lookup)
STATION_WMO = {
    "OUN": "72357", "FWD": "72249", "AMA": "72363", "DDC": "72451",
    "TOP": "72456", "SGF": "72440", "LZK": "72340", "SHV": "72248",
    "BMX": "72230", "JAN": "72235", "LCH": "72240", "BNA": "72327",
    "ILN": "72426", "DTX": "72632", "DVN": "74455", "OAX": "72558",
    "LBF": "72562", "ABR": "72659", "UNR": "72662", "BIS": "72764",
    "GJT": "72476", "DNR": "72469", "ABQ": "72365", "EPZ": "72270",
    "TFX": "72776", "SLC": "72572", "BOI": "72681", "MFR": "72597",
    "OTX": "72786", "UIL": "72797", "REV": "72489", "VBG": "72393",
    "NKX": "72293", "TUS": "72274", "FGZ": "72376", "IAD": "72403",
    "WAL": "72402", "MHX": "72305", "GSO": "72317", "CHS": "72208",
    "JAX": "72206", "TBW": "72210", "MFL": "72202", "TLH": "72214",
    "BUF": "72528", "ALB": "72518", "OKX": "72501", "GYX": "74389",
    "CHH": "74494", "CAR": "72712", "PIT": "72520", "RNK": "72318",
    "MPX": "72649", "GRB": "72645", "ILX": "74560", "APX": "72634",
    "INL": "72747", "CRP": "72251", "MAF": "72265", "DRT": "72261",
    "BRO": "72250", "KEY": "72201", "RIW": "72672", "LMN": "74646",
    "FFC": "72215", "LIX": "72233", "VEF": "72388", "OAK": "72493",
    "GGW": "72768", "ANC": "70273", "FAI": "70261", "LIH": "91165",
    "ITO": "91285",
}


# Available BUFKIT models (for --model flag)
BUFKIT_MODELS = {
    "rap":     "Rapid Refresh (RAP) - hourly, CONUS 13 km",
    "hrrr":    "High-Res Rapid Refresh (HRRR) - hourly, CONUS 3 km",
    "nam":     "North American Mesoscale (NAM) - hourly, 12 km",
    "namnest": "NAM Nest - hourly, 3 km CONUS",
    "gfs":     "Global Forecast System (GFS) - 3-hourly, global",
    "sref":    "Short-Range Ensemble Forecast (SREF) - 3-hourly",
}


# ─── PSU BUFKIT FEED ────────────────────────────────────────────────
PSU_MODELS = {
    "rap":     "Rapid Refresh (RAP)",
    "nam":     "NAM 12 km",
    "namnest": "NAM Nest 3 km",
    "gfs":     "GFS global",
    "hrrr":    "HRRR 3 km",
    "nam4km":  "NAM 4-km CONUS",
    "hiresw":  "HiResW NMMB / ARW",
    "sref":    "SREF ensemble mean",
}



DATA_SOURCES = {
    "obs":    "Observed radiosonde (IEM / UWyo)",
    "rap":    "RAP model analysis - any lat/lon, CONUS (NCEI THREDDS)",
    "bufkit": "BUFKIT forecast models - station-based (Iowa State)",
    "psu":    "PSU BUFKIT feed - latest run, station-based (Penn State)",
    "acars":  "ACARS/AMDAR aircraft obs - airport-based (IEM)",
}


TORNADO_SCAN_STATIONS = [
    "OUN", "FWD", "DDC", "AMA", "TOP", "LMN", "OAX", "SGF", "LZK",
    "SHV", "JAN", "BMX", "ILX", "DVN", "MPX", "BNA", "LBF", "MAF",
    "CRP", "DRT", "BRO", "ABQ", "DNR", "GRB", "APX", "ILN", "RNK",
    "IAD", "MHX", "TLH", "TBW", "JAX",
]

