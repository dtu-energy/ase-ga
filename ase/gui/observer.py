import warnings
import weakref


class Observers:
    def __init__(self):
        self.observer_weakrefs = []

    def register(self, observer):
        if hasattr(observer, '__self__'):  # observer is an instance method
            # Since bound methods are shortlived we need to store the instance
            # it is bound to and use getattr() later:
            obj = observer.__self__
            name = observer.__name__
        else:
            obj = observer
            name = None
        self.observer_weakrefs.append((weakref.ref(obj), name))

    def notify(self):
        # We should probably add an event class to these callbacks.
        weakrefs_still_alive = []
        for weak_ref, name in self.observer_weakrefs:
            observer = weak_ref()
            if observer is not None:
                weakrefs_still_alive.append((weak_ref, name))
                if name is not None:
                    # If the observer is an instance method we stored
                    # self, for garbage collection reasons, and now need to
                    # get the actual method:
                    observer = getattr(observer, name)

                try:
                    observer()
                except Exception as ex:
                    import traceback

                    tb = ''.join(traceback.format_exception(ex))
                    warnings.warn(
                        f'Suppressed exception in observer {observer}: {tb}'
                    )
                    continue

        self.observer_weakrefs = weakrefs_still_alive
