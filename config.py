from pydantic import BaseModel
from pathlib import Path


class Paths(BaseModel):
    local_personal_directory: str
    project_name: str
    sub_project_name: str
    structures: str
    remote_host_name: str
    remote_personal_directory: str


class MetaParameters(BaseModel):
    system_name: str
    type_of_run: str
    parameter_set: str
    num_geometry_steps: int
    geometry_guess_quality: str
    spin_polarized: bool
    continuation: bool
    overwrite_parameters: dict


class POTCARParameters(BaseModel):
    setups: str
    xc: str
    pseudo_potential: str


class KPOINTParameters(BaseModel):
    explicit_kpoints: list[int]


class JobParameters(BaseModel):
    computer_name: str
    job_minutes: int
    num_job_nodes: int
    executable: str
    type_of_run: str  # TODO: This is a duplicate of the type_of_run in MetaParameters
    system_name: str  # TODO: This is a duplicate of the system_name in MetaParameters
    job_memory: int
    job_cores: int
    job_partition: str
    job_qos: str
    job_mail_address: str


class Config(BaseModel):
    submit: bool
    paths: Paths
    meta_parameters: MetaParameters
    potcar_parameters: POTCARParameters
    kpoint_parameters: KPOINTParameters
    dimer_parameters: dict
    job_parameters: JobParameters
