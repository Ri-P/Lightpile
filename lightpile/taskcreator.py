# -*- coding: utf-8 -*-
"""
Based on a dictionary containing only basic data types the environment (
materials, emitters, stacks) objects are initialized, tasks identified and
initialized.
"""
from __future__ import print_function
from __future__ import unicode_literals

from tasks import DipolestudyTask
from lightpile.resources import get_resource_path
from lightpile.material import Material
from lightpile.emitters import OrientedEmitter
from lightpile.stacks import PlanarStack, SingleEmissionplaneStack, PlanarWG
from lightpile.stacks import StackItemError

__author__ = "Richard Pfeifer"


class Taskcreator(object):
    """
    Reads 'user_conf_dict', builds environment and calls Taskinitializers.

    Environment: - materials
                 - emitters
                 - stack
        Environment sections in 'user_conf_dict' are interpreted here and
        environment objects are build.

    All other sections are relayed to the Tasks.
    """
    registered_tasks = {'dipolestudy':  DipolestudyTask}

    def __init__(self, user_conf_dict):
        self.user_conf_dict = user_conf_dict.copy()

    def create_tasks(self):
        """
        Return list of initialized tasks ready to run.
        """
        tasknames = self._identify_tasknames()
        if not len(tasknames) > 0:
            print("No tasks given.")
            return []
        else:
            self._build_environment()
            return [self._create_task(taskname) for taskname in tasknames]

    def _identify_tasknames(self):
        return [sectionname for sectionname in self.user_conf_dict.keys()
                if sectionname in self.registered_tasks]

    def _build_environment(self):
        self._build_materials()
        self._build_emitters()
        self._build_stack()

    def _create_task(self, taskname):
        """
        Find task-class from taskname and initialize task-object
        with environment objects and userdata.
        """
        Task = self.registered_tasks[taskname]
        return Task(materials=self.materials, emitters=self.emitters,
                    stack=self.stack, **self.user_conf_dict)

    def _build_materials(self):
        """
        Build 'self.materials' from 'materials' section in user_conf_dict.
        """
        self.materials = {}
        for material_data in self.user_conf_dict.pop("materials"):
            if material_data["name"] in self.materials:
                raise ValueError(
                    "Materialname {0} is already known."
                    "Every materialname has to be unique. "
                    "Did you name two materials the same?"
                    .format(material_data["name"]))
            if isinstance(material_data["nk"], str):
                material_data["nk"] = get_resource_path(material_data["nk"])
            m = Material(material_data["nk"])
            self.materials.update({material_data["name"]: m})

    def _build_emitters(self):
        """
        Build 'self.emitters' from 'emitters' section in user_conf_dict.

        Only emitterdescriptions matching 'OrientedEmitter' are accepted.
        """
        self.emitters = {}
        emitter_descriptions = self.user_conf_dict.pop("emitters")
        if len(emitter_descriptions) != 1:
            raise ValueError("Error. Expected exactly one emitter.")
        for emitter_description in emitter_descriptions:
            for key in emitter_description.keys():
                if key not in ["name", "orientation", "iqe"]:
                    raise ValueError("Unknown emitter property {0}.".format(
                        key))
            if emitter_description["name"] in self.emitters:
                raise ValueError(
                    "Emittername {0} is already known."
                    "Every emittername has to be unique. "
                    "Did you name two emitters the same?"
                    .format(emitter_description["name"]))
            e = OrientedEmitter(
                orientation=emitter_description["orientation"],
                iqe=emitter_description["iqe"]
            )
            self.emitters.update({emitter_description["name"]: e})

    def _build_stack(self):
        """
        Read 'stack' section and build 'self.stack'.
        """
        stacked_objects = []
        has_emitter = False
        for stacked_token in self.user_conf_dict.pop("stack"):
            if "planarmaterial" in stacked_token:
                try:
                    mat = self.materials[stacked_token["planarmaterial"]]
                except KeyError:
                    raise StackItemError(
                        "Error: Material '{0}' in stack "
                        "was not declared in section 'Materials'."
                        .format(stacked_token["planarmaterial"]))
                else:
                    stacked_objects.append(
                        PlanarWG(mat, stacked_token["thickness"]))
            elif "emittername" in stacked_token:
                name = stacked_token["emittername"]
                try:
                    stacked_objects.append(self.emitters[name])
                except KeyError:
                    raise Exception(
                        "Emitter {0} not defined in 'emitters' "
                        "section.".format(name))
                else:
                    has_emitter = True
            else:
                raise Exception(
                    "Unknown stack token '{0}'"
                    .format(stacked_token.__repr__()))
        if has_emitter:
            self.stack = SingleEmissionplaneStack(stacked_objects)
        else:
            self.stack = PlanarStack(stacked_objects)
