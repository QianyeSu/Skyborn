/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida
          Oliver Heidmann

*/
#ifndef TEXT_H
#define TEXT_H

#include <string>
#include <sstream>
#include <type_traits>

enum ColorEnabled
{
  No,
  All,
  Auto
};

void mpmo_color_set(int p_colorBehaviour);

enum TextMode
{
  MODELESS = 0,
  BRIGHT = 1,
  DIM = 2,
  UNDERLINE = 4,
  BLINK = 5,
  REVERSE = 7,
  HIDDEN = 8
};

// Used for text color and background
enum TextColor
{
  NO_COLOR = 0,
  BLACK = 30,
  RED = 31,
  GREEN = 32,
  YELLOW = 33,
  BLUE = 34,
  MAGENTA = 35,
  CYAN = 36,
  WHITE = 37
};

bool color_enabled();

int mpmo_get_color_mode();

void set_text_color(std::FILE *fp, TextMode attr, TextColor fg);
void set_text_color(std::FILE *fp, TextMode attr);
void set_text_color(std::FILE *fp, TextColor fg);
void reset_text_color(std::FILE *fp);

// Templates for converting any given input to a string
template <class T>
typename std::enable_if<std::is_fundamental<T>::value, std::string>::type
stringify(T &t)
{
  return std::to_string(t);
}

template <class T>
typename std::enable_if<!std::is_fundamental<T>::value, std::string>::type
stringify(T &t)
{
  return std::string(t);
}

// public color functions

/* function that returns a String containing a ANSI escape code for color, background color and specific modes.
 * for reference look up 'ANSI escape code'
 */
static inline std::string
text_color(TextColor foreground = NO_COLOR, TextMode mode = MODELESS, TextColor background = NO_COLOR)
{
  bool tty = true;
  std::stringstream s;

  if (!color_enabled()) return "";
  s << "\033[";
  if (!tty)
  {
    s << "m";
    return s.str();
  }
  if (!foreground && !background)
  {
    s << "0";  // reset colors if no params
  }
  if (foreground)
  {
    s << foreground;
    if (background) s << ";";
  }
  if (background)
  {
    s << 10 + background;
    if (mode) s << ";";
  }
  else if (mode) { s << ";"; }
  if (mode) { s << mode; }
  s << "m";
  return s.str();
}

/********************************************************************************/

// Templates for colorizing any variable or string

template <class T>
std::string
Red(T p_str)
{
  return text_color(RED) + stringify(p_str) + text_color();
}

template <class T>
std::string
Cyan(T p_str)
{
  return text_color(CYAN) + stringify(p_str) + text_color();
}

template <class T>
std::string
Blue(T p_str)
{
  return text_color(BLUE) + stringify(p_str) + text_color();
}

template <class T>
std::string
Green(T p_str)
{
  return text_color(GREEN) + stringify(p_str) + text_color();
}

template <class T>
std::string
Yellow(T p_str)
{
  return text_color(YELLOW) + stringify(p_str) + text_color();
}

template <class T>
std::string
Black(T p_str)
{
  return text_color(BLACK) + stringify(p_str) + text_color();
}

/********************************************************************************/

#endif
