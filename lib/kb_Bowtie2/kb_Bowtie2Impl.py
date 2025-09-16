# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
from pprint import pprint
from kb_Bowtie2.util.Bowtie2IndexBuilder import Bowtie2IndexBuilder
from kb_Bowtie2.util.Bowtie2Aligner import Bowtie2Aligner
from kb_Bowtie2.util.Bowtie2Runner import Bowtie2Runner
#END_HEADER


class kb_Bowtie2:
    '''
    Module Name:
    kb_Bowtie2

    Module Description:
    A KBase module: kb_Bowtie2
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.1.0"
    GIT_URL = "git@github.com:kbaseapps/kb_Bowtie2.git"
    GIT_COMMIT_HASH = "fc2f147ea45b59fa54ab2663b06b40043613ffc5"

    #BEGIN_CLASS_HEADER
    @staticmethod
    def _truthy(val):
        """Return True if val represents an enabled checkbox/bool-like value."""
        if isinstance(val, bool):
            return val
        if isinstance(val, (int, float)):
            return val != 0
        if isinstance(val, str):
            return val.strip().lower() in ("1", "true", "t", "yes", "y", "on")
        return False

    @staticmethod
    def _as_int_or_none(val):
        if val is None:
            return None
        try:
            return int(val)
        except Exception:
            return None

    @classmethod
    def _normalize_params(cls, params, for_index=False):
        """
        Lightly sanitize/normalize user params so downstream util classes
        can consume consistent types/keys.
        - Drops empty sam_opt_config so older code won't emit an empty flag.
        - Normalizes use_sais to a canonical boolean (and also 0/1 int for convenience).
        - Best-effort int coercion for numeric text fields.
        """
        if not isinstance(params, dict):
            return params

        p = dict(params)  # shallow copy

        # --- Bowtie2 ≥2.5 advanced flags ---
        # Keep only non-empty sam_opt_config
        sam_cfg = p.get("sam_opt_config")
        if isinstance(sam_cfg, str):
            sam_cfg = sam_cfg.strip()
        if sam_cfg:
            p["sam_opt_config"] = sam_cfg
        else:
            p.pop("sam_opt_config", None)

        # Normalize use_sais (used by index builder)
        if "use_sais" in p:
            enabled = cls._truthy(p.get("use_sais"))
            p["use_sais"] = enabled
            # some downstream code prefers 0/1; provide both
            p["use_sais_int"] = 1 if enabled else 0

        # Coerce common numeric text fields to ints when possible (harmless if unused)
        for k in ("trim5", "trim3", "np", "minins", "maxins"):
            if k in p:
                iv = cls._as_int_or_none(p.get(k))
                if iv is not None:
                    p[k] = iv

        # Alignment type / quality score: strip whitespace for safety
        for k in ("alignment_type", "quality_score", "preset_options", "orientation", "condition_label"):
            if k in p and isinstance(p[k], str):
                p[k] = p[k].strip()

        # Optionally mark this call's intent
        if for_index:
            p["_intent"] = "build_index"
        else:
            p["_intent"] = "align"

        return p
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch_dir = os.path.abspath(config['scratch'])
        self.workspace_url = config['workspace-url']
        self.srv_wiz_url = config['srv-wiz-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        #END_CONSTRUCTOR
        pass


    def align_reads_to_assembly_app(self, ctx, params):
        """
        :param params: instance of type "AlignReadsParams" ...
        :returns: instance of type "AlignReadsResult" ...
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN align_reads_to_assembly_app
        print('Running align_reads_to_assembly_app() with raw params=')
        pprint(params)

        # Normalize & pass through the new advanced flags
        nparams = self._normalize_params(params, for_index=False)

        # Log the normalized subset of interest for debug
        debug_view = {
            k: nparams.get(k)
            for k in (
                "alignment_type", "quality_score", "preset_options", "trim5", "trim3",
                "np", "minins", "maxins", "orientation",
                # new:
                "sam_opt_config"
            )
        }
        print('Normalized params for align (subset):')
        pprint(debug_view)

        bowtie2_aligner = Bowtie2Aligner(self.scratch_dir, self.workspace_url,
                                         self.callback_url, self.srv_wiz_url,
                                         ctx.provenance())
        result = bowtie2_aligner.align(nparams)
        #END align_reads_to_assembly_app

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method align_reads_to_assembly_app return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def align_one_reads_to_assembly(self, ctx):
        """
        aligns a single reads object to produce
        """
        # ctx is the context object
        #BEGIN align_one_reads_to_assembly
        #END align_one_reads_to_assembly
        pass

    def get_bowtie2_index(self, ctx, params):
        """
        :param params: instance of type "GetBowtie2Index" ...
        :returns: instance of type "GetBowtie2IndexResult" ...
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN get_bowtie2_index
        print('Running get_bowtie2_index() with raw params=')
        pprint(params)

        nparams = self._normalize_params(params, for_index=True)

        # Show whether SAIS is enabled after normalization
        print('Index build flags (subset): {"use_sais": %s, "use_sais_int": %s}' %
              (nparams.get("use_sais"), nparams.get("use_sais_int")))

        bowtie2IndexBuilder = Bowtie2IndexBuilder(self.scratch_dir, self.workspace_url,
                                                  self.callback_url, self.srv_wiz_url,
                                                  ctx.provenance())
        result = bowtie2IndexBuilder.get_index(nparams)
        #END get_bowtie2_index

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method get_bowtie2_index return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def run_bowtie2_cli(self, ctx, params):
        """
        general purpose local function for running tools in the bowtie2 suite
        :param params: instance of type "RunBowtie2CLIParams" ...
        """
        # ctx is the context object
        #BEGIN run_bowtie2_cli
        print('Running run_bowtie2_cli() with params=')
        pprint(params)

        # accept either 'command' or 'command_name'
        command = params.get('command', params.get('command_name'))
        options = params.get('options')

        if not command:
            raise ValueError('required parameter "command" (or "command_name") was missing.')
        if options is None:
            raise ValueError('required parameter field "options" was missing.')

        bowtie2 = Bowtie2Runner(self.scratch_dir)
        bowtie2.run(command, options)
        #END run_bowtie2_cli
        pass

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
