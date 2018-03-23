NAME = nlt

CFILES = nlt.c adlift.c

LIBS = -lm

CFLAGS = -Wall -Wextra -Wfloat-equal -pedantic -std=c99

DFLAGS = -g

OFLAGS = -O3

$(NAME) : $(CFILES)
	#gcc $(CFILES) -o $(NAME) $(OFLAGS) $(CFLAGS) $(LIBS)
	gcc $(CFILES) -o nlt -I/usr/share/R/include -DNDEBUG -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-rAT5Oi/r-base-3.3.1=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g $(CFLAGS) $(LIBS)

run: $(NAME)
	./$(NAME)

clean:
	rm $(NAME)

profile:
	gcc $(CFILES) -o $(NAME) $(OFLAGS) $(CFLAGS) $(LIBS) -pg

profileout:
	gprof $(NAME) gmon.out > analysis.txt

debug:
	gcc $(CFILES) -o $(NAME) $(OFLAGS) $(DFLAGS) $(CFLAGS) $(NO_WARN_UNUSED) $(LIBS)

rungdb:
	gdb $(NAME)

memcheck:
	valgrind --tool=memcheck --leak-check=yes ./$(NAME)

test: $(NAME)
	./$(NAME) best made
