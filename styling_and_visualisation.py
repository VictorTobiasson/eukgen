#main top level configuration for altair plots, default 
import altair as alt

# color specification
# orange = ['#fff2e6ff',  '#ffe0caff',  '#ffcea9ff',  '#fac093ff',  '#f5ad76ff',  '#eb8f49ff',  '#e4751fff']
# yellow = ['#fffde9ff', '#fef9c4ff', '#fef58bff', '#f8ea6bff', '#f5e53cff', '#f1da09ff', '#f0ce0aff']
# moss = ['#f2ffe8ff', '#e3ffc4ff', '#d1fca5ff', '#b8f775ff', '#a7f258ff', '#8be332ff', '#7ed82fff']
# green = ['#e4ffebff', '#c2ffd2ff', '#a1ffb9ff', '#8dfaa9ff', '#5cef82ff', '#33cf5bff', '#2eb851ff']
# #ice = ['#e6f8faff', '#c5f4faff', '#a7f1faff', '#84e9f5ff', '#74d8e3ff', '#2bbecfff', '#09a1b2ff']
# #blue = ['#eff3ffff', '#dee7ffff', '#c4d7ffff', '#a7beffff', '#8ba9ffff', '#5c7fe1ff', '#3f69e0ff']
# #purple = ['#f7f0ffff', '#eddbffff', '#e0c3ffff', '#d2a7ffff', '#ba85f3ff', '#a058e9ff', '#8c3fdeff']
# #pink = ['#fcf0f8ff', '#f7d2eeff', '#f4a8e2ff', '#e984d1ff', '#e160c2ff', '#cf3dacff', '#b63898ff']
# #red = ['#ffebebff', '#fad2d2ff', '#f6b8b8ff', '#f6aaaaff', '#f68f8fff', '#e56969ff', '#ca5353ff']
# #lightness = 6
# spectral9 = [orange[lightness], blue[lightness], red[lightness], green[lightness], pink[lightness], moss[lightness], purple[lightness], yellow[lightness], ice[lightness]] 

warmgrays = ["#cec5c1", "#c0b8b4", "#b3aaa7", "#a59c99", "#98908c", "#8b827f", "#7e7673", "#726866", "#665c5a", "#4c4443"]

light_gray = warmgrays[0]
med_gray = warmgrays[2]
dark_gray = warmgrays[9]
white = '#ffffff'
black = '#000000'

# from seaborn twilight_shifted_r then shifted to every 4th element with an offset of 1
# creats a categorical map from a sequential diverging
# twilight_shifted_r_perm = ['#7f2350', '#d6c2b6', '#6989be', '#501444', '#ca997c', '#89adc5', '#491564', '#bc6b59', '#bccbd1', '#5c359a', '#a54350', '#e2d9e2', '#5f61b4']

colorlib = {'orange': ['#fff2e6ff',  '#ffe0caff',  '#ffcea9ff',  '#fac093ff',  '#f5ad76ff',  '#eb8f49ff',  '#e4751fff'],
        'yellow': ['#fffde9ff', '#fef9c4ff', '#fef58bff', '#f8ea6bff', '#f5e53cff', '#f1da09ff', '#f0ce0aff'],
        'moss': ['#f2ffe8ff', '#e3ffc4ff', '#d1fca5ff', '#b8f775ff', '#a7f258ff', '#8be332ff', '#7ed82fff'],
        'green': ['#e4ffebff', '#c2ffd2ff', '#a1ffb9ff', '#8dfaa9ff', '#5cef82ff', '#33cf5bff', '#2eb851ff'],
        'ice': ['#e6f8faff', '#c5f4faff', '#a7f1faff', '#84e9f5ff', '#74d8e3ff', '#2bbecfff', '#09a1b2ff'],
        'blue': ['#eff3ffff', '#dee7ffff', '#c4d7ffff', '#a7beffff', '#8ba9ffff', '#5c7fe1ff', '#3f69e0ff'],
        'purple': ['#f7f0ffff', '#eddbffff', '#e0c3ffff', '#d2a7ffff', '#ba85f3ff', '#a058e9ff', '#8c3fdeff'],
        'pink': ['#fcf0f8ff', '#f7d2eeff', '#f4a8e2ff', '#e984d1ff', '#e160c2ff', '#cf3dacff', '#b63898ff'],
        'red': ['#ffebebff', '#fad2d2ff', '#f6b8b8ff', '#f6aaaaff', '#f68f8fff', '#e56969ff', '#ca5353ff'],
        'warmgrays': warmgrays, 
        'light_gray': light_gray, 
        'med_gray': med_gray, 
        'dark_gray': dark_gray, 
        'twilight_shifted_r_perm':['#983550', '#cca389', '#c4ced4', '#6276ba', '#b25652', '#d8c7be', '#95b5c7', '#5e51ad', '#c27c63', '#e2d9e2', '#7297c1']}

def default_theme(colorlib=colorlib):
    header_font = 'Nunito'
    body_font = 'NanumGothic'
    body_font_bold = 'NanumGothic'
    base_color = colorlib['twilight_shifted_r_perm'][5]
    
    main_style = {
        "config": {
            "view": {
                "width": 700,
                "height": 525
            },
            "title": {
                "fontSize": 24,
                "font": header_font,
                "fontWeight": 400,
                "anchor": "start", # equivalent of left-aligned
                "color": black
            },
            "axisX": {
                "domain": True,
                "domainColor": colorlib['dark_gray'],
                "domainWidth": 1,
                "grid": True,
                "gridColor": colorlib['light_gray'],
                "gridWidth": 0.5,
                "labelFont": body_font,
                "labelFontSize": 12,
                "labelColor": colorlib['dark_gray'],
                "labelAngle": 0, 
                "labelFlush": False,
                "tickColor": colorlib['dark_gray'],
                "tickSize": 5,
                "titleFont": body_font,
                "titleFontWeight": 100,
                "titleColor": colorlib['dark_gray'],
                "titleFontSize": 14,
                "titlePadding": 10
            },
            "axisY": {
                "domain": True,
                "domainColor": colorlib['dark_gray'],
                "domainWidth": 1,
                "grid": True,
                "gridColor": colorlib['light_gray'],
                "gridWidth": 0.5,
                "labelFont": body_font,
                "labelFontSize": 12,
                "labelColor": colorlib['dark_gray'],
                "labelAngle": 0,
                "tickColor": colorlib['dark_gray'],
                "tickSize": 5,
                "titleFont": body_font,
                "titleFontWeight": 100,
                "titleColor": colorlib['dark_gray'],
                "titleFontSize": 14,
                "titlePadding": 10
            },
            "header": {
                "labelFont": body_font,
                "labelFontSize": 16,
                "titleFont": body_font,
                "titleFontSize": 16
            },
            "range": {
                "category": colorlib['twilight_shifted_r_perm'],
                "diverging": colorlib['twilight_shifted_r_perm'],
                "sequential": colorlib['twilight_shifted_r_perm'],
            },
            "legend": {
                "labelFont": body_font,
                "labelFontSize": 12,
                "symbolSize": 100, # default,
                "titleFont": body_font,
                "titleFontSize": 12,
            },
            ### MARKS CONFIGURATIONS ###
            "area": {
               "fill": base_color,
            },
            "circle": {
               "fill": base_color,
               "size": 40,
               "opacity": 1,
               "strokeWidth": 0.2,
               "stroke": colorlib['dark_gray'],
               "scheme": 'category20b',
            },
            "line": {
               "color": base_color,
               "stroke": base_color,
               "strokeWidth": 1.5,
            },
            "trail": {
               "color": base_color,
               "stroke": base_color,
               "strokeWidth": 0,
               "size": 1,
            },
            "path": {
               "stroke": base_color,
               "strokeWidth": 0.5,
            },
            "point": {
               "color": base_color,
               "size": 40,
               "strokeWidth": 0.5,
            },
            "text": {
               "font": body_font,
               "color": base_color,
               "fontSize": 11,
               "align": "right",
               "size": 14,            
            }, 
            "rule": {
               "color": colorlib['med_gray'],
            }, 
            "bar": {
                #"size": 10,
                "binSpacing": 1,
                "continuousBandSize": 10,
#                 "discreteBandSize": 10,
                "fill": base_color,
                "stroke": False,
            },
            "tick": {
                "color": base_color
            }
        }   
    }

    return main_style

def set_alt_theme():
    alt.themes.register("default_theme", default_theme)
    alt.themes.enable("default_theme")

if __name__ == '__main__':
    set_alt_theme()

