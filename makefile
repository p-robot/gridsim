CC=g++
CFLAGS=-c -std=gnu++11 -O2
INCLUDE=-I./include/
OBJDIR=obj/

OBJLIST   = Cell.o common_functions.o Config.o Grid.o Local_spread.o Node.o Point.o Shipment_handler.o main.o

OBJECTS   = $(addprefix $(OBJDIR), $(OBJLIST) )

all:gridsim

gridsim: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

$(OBJECTS): ./$(OBJDIR)%.o: src/%.cpp
	$(CC) $(CFLAGS) $? -o $@ $(INCLUDE)

clean:
	rm obj/*.o
