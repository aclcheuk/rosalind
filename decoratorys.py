def announce(f):
    def wrapper():
        print("About to start function...")
        f()
        print("Completed function")
    return wrapper

@announce
def hello():
    print("This is the function step")

hello()