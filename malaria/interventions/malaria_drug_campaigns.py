from malaria.interventions.malaria_drugs import drug_configs_from_code
from malaria.interventions.malaria_diagnostic import add_diagnostic_survey  # , add_triggered_survey
from copy import deepcopy, copy
import random


def add_drug_campaign(cb, campaign_type, drug_code, start_days, coverage=1.0, repetitions=3, interval=60,
                      diagnostic_threshold=40, fmda_radius='hh', node_selection_type='DISTANCE_ONLY',
                      trigger_coverage=1.0, snowballs=0, treatment_delay=0, triggered_campaign_delay=0, nodes=[],
                      target_group='Everyone', dosing='', drug_ineligibility_duration=0,
                      node_property_restrictions=[], ind_property_restrictions=[], trigger_condition_list=[],
                      listening_duration=-1):
    """
    Add a drug campaign defined by the parameters to the config builder.
    Note: When using "trigger_condition_list", the first entry of "start_days" is the day that is used to start
    listening for the trigger(s), if listening_duration is specified, the rest of the days are used to start others d;
     "triggered_campaign_delay" indicates the delay period between receiving the trigger and the first campaign/intervention.
     "listening_duration" is the duration for which the listen for the trigger, -1 indicates "indefinitely/forever".

    :param cb: The :py:class:`DTKConfigBuilder <dtk.utils.core.DTKConfigBuilder>` that will receive the drug
    intervention.
    :param campaign_type: type of drug campaign (MDA, MSAT, SMC, fMDA, rfMSAT, or rfMDA)
    # fMDA, rfMDA, rfMSAT are appropriate for household-type simulations only.
    :param drug_code: The drug code of the drug regimen (AL, DP, etc; allowable types defined in malaria_drugs.py)
    :param start_days: List of start days where the drug regimen will be distributed
    :param coverage: Demographic coverage of the distribution (fraction of people at home during campaign)
    :param repetitions: Number repetitions
    :param interval: Timesteps between the repetitions. For RCD (rfMDA, rfMSAT), interval indicates the duration of RCD
    :param diagnostic_threshold: diagnostic sensitivity for diagnostic-dependent campaigns (MSAT, fMDA, rfMSAT)
    :param fmda_radius: radius of focal response upon finding infection, in km or 'hh' for within-household only
    :param node_selection_type: restriction on broadcasting focal response trigger.
    # DISTANCE_ONLY: It will send the event to nodes that are within a given distance.
    # MIGRATION_NODES_ONLY: It will only send the event to nodes that the individual can migrate to.
    # DISTANCE_AND_MIGRATION: It will only send the even to migratable nodes that are within a given distance.
    # Migrateable nodes = Local and Regional.
    :param trigger_coverage: for RCD only. Fraction of trigger events that will trigger an RCD. "coverage" param sets
    the fraction of individuals reached during RCD response.
    :param snowballs: number of snowball levels in reactive response.
    :param treatment_delay: for MSAT and fMDA, length of time between administering diagnostic and giving drugs; for RCD, length
    of time between treating index case and triggering RCD response.
    :param triggered_campaign_delay: when using trigger_condition_list, the delay between receiving the trigger event and
    running the triggered campaign/intervention
    :param nodes: list of node IDs; if empty, defaults to all nodes
    :param target_group: dict of {'agemin' : x, 'agemax' : y} to target MDA, SMC, MSAT, fMDA to individuals between
    x and y years of age.
    For fMDA, only the testing portion is targeted; everyone is eligible for treating.
    :param dosing: change drug dosing style. Default is FullTreatmentCourse.
    :param drug_ineligibility_duration: if this param is > 0, use IndividualProperties to prevent people from receiving
    drugs too frequently. Demographics file will need to define the IP DrugStatus with possible values None and
    RecentDrug. Individuals with status RecentDrug will not receive drugs during drug campaigns, though they are still
    eligible for receiving diagnostics (in MSAT, etc). Individuals who receive drugs during campaigns will have their
    DrugStatus changed to RecentDrug for drug_ineligibility_duration days.
    :param node_property_restrictions: used with NodePropertyRestrictions.
    :param trigger_condition_list: When trigger_string is set, the first entry of start_days is the day that is used to start
    listening for the trigger(s), the campaign happens when the trigger(s) is received.
    :param listening_duration: is the duration for which the listen for the trigger, -1 indicates "indefinitely/forever"
    Format: list of dicts: [{ "NodeProperty1" : "PropertyValue1" }, {'NodeProperty2': "PropertyValue2"}, ...]
    """

    expire_recent_drugs = {}
    if drug_ineligibility_duration > 0:
        expire_recent_drugs = {"class": "PropertyValueChanger",
                               "Target_Property_Key": "DrugStatus",
                               "Target_Property_Value": "RecentDrug",
                               "Daily_Probability": 1.0,
                               "Maximum_Duration": 0,
                               'Revert': drug_ineligibility_duration
                               }

    # set up intervention drug block
    drug_configs = drug_configs_from_code(cb, drug_code)

    # update drug dosing if requested
    if dosing != '':
        for i in range(len(drug_configs)):
            drug_configs[i]['Dosing_Type'] = dosing

    # set up events to broadcast when receiving campaign drug
    receiving_drugs_event = {
        "class": "BroadcastEvent",
        "Broadcast_Event": "Received_Campaign_Drugs"
    }
    if 'Vehicle' in drug_code:  # if distributing Vehicle drug
        receiving_drugs_event["Broadcast_Event"] = "Received_Vehicle"
    if campaign_type[0] == 'r':  # if reactive campaign
        receiving_drugs_event['Broadcast_Event'] = 'Received_RCD_Drugs'

    # set up node config block
    node_cfg = {"class": "NodeSetAll"}
    if nodes:
        node_cfg = {'Node_List': nodes, "class": "NodeSetNodeList"}

    # set up drug campaign
    if campaign_type == 'MDA' or campaign_type == 'SMC':  # standard drug campaign: MDA, no event triggering
        add_MDA(cb, start_days, coverage, drug_configs, receiving_drugs_event, repetitions, interval, node_cfg,
                expire_recent_drugs, node_property_restrictions, ind_property_restrictions, target_group,
                trigger_condition_list, listening_duration, triggered_campaign_delay)

    elif campaign_type == 'MSAT':  # standard drug campaign: MSAT, no event triggering
        add_MSAT(cb, start_days, coverage, drug_configs, receiving_drugs_event, repetitions, interval,
                 treatment_delay, diagnostic_threshold, node_cfg, expire_recent_drugs, node_property_restrictions,
                 ind_property_restrictions, target_group, trigger_condition_list, triggered_campaign_delay,
                 listening_duration)

    elif campaign_type == 'fMDA':
        add_fMDA(cb, start_days, trigger_coverage, coverage, drug_configs, receiving_drugs_event, repetitions, interval,
                 treatment_delay, diagnostic_threshold, fmda_radius, node_selection_type, node_cfg, expire_recent_drugs,
                 node_property_restrictions, ind_property_restrictions, target_group, trigger_condition_list,
                 listening_duration,  triggered_campaign_delay)

    elif campaign_type == 'rfMSAT':
        add_rfMSAT(cb, start_days[0], coverage, drug_configs, receiving_drugs_event, interval, treatment_delay,
                   trigger_coverage, diagnostic_threshold, fmda_radius, node_selection_type, snowballs, node_cfg,
                   expire_recent_drugs, node_property_restrictions, ind_property_restrictions, trigger_condition_list,
                   listening_duration, triggered_campaign_delay)

    elif campaign_type == 'rfMDA':
        add_rfMDA(cb, start_days[0], coverage, drug_configs, receiving_drugs_event, interval, treatment_delay,
                  trigger_coverage, fmda_radius, node_selection_type, node_cfg, expire_recent_drugs,
                  node_property_restrictions, ind_property_restrictions, trigger_condition_list, listening_duration,
                  triggered_campaign_delay)

    elif campaign_type in ['borderscreen']:
        # PROBABLY DOES NOT WORK --- WANT ONLY PEOPLE COMING IN FROM WORK NODE,
        # NOT DEPARTING OR RETURNING FROM LOCAL NODES
        add_borderscreen(cb, start_days[0], coverage, drug_configs, receiving_drugs_event, interval, treatment_delay,
                         trigger_coverage, diagnostic_threshold, fmda_radius, node_selection_type, snowballs, node_cfg,
                         expire_recent_drugs, trigger_condition_list, triggered_campaign_delay, listening_duration)

    else:
        pass

    return {'drug_campaign.type': campaign_type,
            'drug_campaign.drugs': drug_code,
            'drug_campaign.trigger_coverage': trigger_coverage,
            'drug_campaign.coverage': coverage
            }


def add_MDA(cb, start_days, coverage, drug_configs, receiving_drugs_event, repetitions, interval,
            nodes, expire_recent_drugs, node_property_restrictions, ind_property_restrictions, target_group,
            trigger_condition_list=[], listening_duration=-1, triggered_campaign_delay=0):
    if trigger_condition_list:
        sorta_unique_id = random.randrange(10000)
        for repetition in range(0, repetitions):
            # set up trigger event that triggers every person to send out the the trigger for the repetition's listeners
            if repetition == 0 and repetitions > 1:
                trigger_event = {
                    "class": "CampaignEvent",
                    "Start_Day": start_days[0],
                    "Nodeset_Config": nodes,
                    "Event_Coordinator_Config": {
                        "class": "StandardInterventionDistributionEventCoordinator",
                        "Intervention_Config": {
                            "class": "NodeLevelHealthTriggeredIV",
                            "Trigger_Condition_List": trigger_condition_list,
                            "Duration": listening_duration,
                            "Target_Residents_Only": 1,
                            "Actual_IndividualIntervention_Config": {
                                "class": "MultiInterventionDistributor",
                                "Intervention_List": []
                            }
                        }
                    }
                }
                for x in range(1, repetitions):
                    delayed_later_event_trigger = {
                        "class": "DelayedIntervention",
                        "Delay_Distribution": "FIXED_DURATION",
                        "Delay_Period": triggered_campaign_delay + interval * x,
                        "Actual_IndividualIntervention_Configs":
                            [
                                {
                                    "class": "BroadcastEvent",
                                    "Broadcast_Event": str(sorta_unique_id) + "_" + str(x)
                                }
                            ]
                    }
                    trigger_event['Event_Coordinator_Config']['Intervention_Config'][
                        "Actual_IndividualIntervention_Config"]["Intervention_List"].append(delayed_later_event_trigger)

                cb.add_event(trigger_event)

            drug_event = {
                "class": "CampaignEvent",
                "Start_Day": start_days[0],
                "Nodeset_Config": nodes,
                "Event_Coordinator_Config":
                    {
                        "class": "StandardInterventionDistributionEventCoordinator",
                        "Intervention_Config":
                            {
                                "class": "NodeLevelHealthTriggeredIV",
                                "Blackout_On_First_Occurrence": 1,
                                "Trigger_Condition_List": trigger_condition_list,
                                "Demographic_Coverage": coverage,
                                "Duration": listening_duration,
                                "Target_Residents_Only": 1,
                                "Actual_IndividualIntervention_Config":
                                    {
                                        "class": "DelayedIntervention",
                                        "Delay_Distribution": "FIXED_DURATION",
                                        "Delay_Period": triggered_campaign_delay,
                                        "Actual_IndividualIntervention_Configs":
                                            [
                                                {
                                                    "class": "MultiInterventionDistributor",
                                                    "Intervention_List": drug_configs + [receiving_drugs_event]
                                                }
                                            ]
                                    }
                            }
                    }
            }

            if repetitions > 1:
                if repetition > 0:
                    drug_event['Event_Coordinator_Config']['Intervention_Config']["Trigger_Condition_List"] = [
                        str(sorta_unique_id) + "_" + str(repetition)]

            if ind_property_restrictions:
                drug_event['Event_Coordinator_Config']['Intervention_Config'][
                    "Property_Restrictions_Within_Node"] = ind_property_restrictions

            if expire_recent_drugs:
                drugstatus = {"DrugStatus": "None"}
                if ind_property_restrictions:
                    drug_event['Event_Coordinator_Config']['Intervention_Config'][
                        "Property_Restrictions_Within_Node"] = [
                        dict(drugstatus.items() + x.items()) for x in ind_property_restrictions]
                else:
                    drug_event['Event_Coordinator_Config']['Intervention_Config'][
                        "Property_Restrictions_Within_Node"] = [drugstatus]
                drug_event['Event_Coordinator_Config']['Intervention_Config']["Actual_IndividualIntervention_Config"][
                    "Actual_IndividualIntervention_Configs"][0]["Intervention_List"].append(expire_recent_drugs)

            if target_group != 'Everyone':
                drug_event['Event_Coordinator_Config']['Intervention_Config'].update({
                    "Target_Demographic": "ExplicitAgeRanges",  # Otherwise default is Everyone
                    "Target_Age_Min": target_group['agemin'],
                    "Target_Age_Max": target_group['agemax']
                })
            else:
                drug_event['Event_Coordinator_Config']['Intervention_Config'][
                    'Actual_IndividualIntervention_Config'].update({"Target_Demographic": "Everyone"})

            if node_property_restrictions:
                drug_event['Event_Coordinator_Config']['Intervention_Config'][
                    'Node_Property_Restrictions'] = node_property_restrictions

            cb.add_event(drug_event)

    else:
        for start_day in start_days:
            drug_event = {
                "class": "CampaignEvent",
                "Start_Day": start_day,
                "Event_Coordinator_Config": {
                    "class": "StandardInterventionDistributionEventCoordinator",
                    "Target_Demographic": "Everyone",
                    "Demographic_Coverage": coverage,
                    "Intervention_Config": {
                        "class": "MultiInterventionDistributor",
                        "Intervention_List": drug_configs + [receiving_drugs_event]
                    },
                    "Number_Repetitions": repetitions,
                    "Timesteps_Between_Repetitions": interval
                },
                "Nodeset_Config": nodes
            }

            drug_event['Event_Coordinator_Config']['Intervention_Config'] = {
                "class": "MultiInterventionDistributor",
                "Intervention_List": drug_configs + [receiving_drugs_event]
            }

            if ind_property_restrictions:
                drug_event['Event_Coordinator_Config']["Property_Restrictions_Within_Node"] = ind_property_restrictions

            if expire_recent_drugs:
                drugstatus = {"DrugStatus": "None"}
                if ind_property_restrictions:
                    drug_event['Event_Coordinator_Config']["Property_Restrictions_Within_Node"] = [
                        dict(drugstatus.items() + x.items()) for x in ind_property_restrictions]
                else:
                    drug_event['Event_Coordinator_Config']["Property_Restrictions_Within_Node"] = [drugstatus]
                drug_event['Event_Coordinator_Config']["Intervention_Config"]["Intervention_List"].append(
                    expire_recent_drugs)

            if target_group != 'Everyone':
                drug_event['Event_Coordinator_Config'].update({
                    "Target_Demographic": "ExplicitAgeRanges",  # Otherwise default is Everyone
                    "Target_Age_Min": target_group['agemin'],
                    "Target_Age_Max": target_group['agemax']
                })
            if node_property_restrictions:
                drug_event['Event_Coordinator_Config']['Node_Property_Restrictions'] = node_property_restrictions
            cb.add_event(drug_event)


def add_MSAT(cb, start_days, coverage, drug_configs, receiving_drugs_event, repetitions, interval,
             treatment_delay, diagnostic_threshold, nodes, expire_recent_drugs, node_property_restrictions,
             ind_property_restrictions, target_group, trigger_condition_list, triggered_campaign_delay, listening_duration):
    event_config = drug_configs + [receiving_drugs_event]
    IP_restrictions = []
    if expire_recent_drugs:
        event_config.append(expire_recent_drugs)
        IP_restrictions = [{"DrugStatus": "None"}]

    if treatment_delay == 0:
        msat_cfg = event_config
    else:
        msat_cfg = [{"class": "DelayedIntervention",
                     "Delay_Distribution": "FIXED_DURATION",
                     "Delay_Period": treatment_delay,
                     "Actual_IndividualIntervention_Configs": event_config
                     }]

    if trigger_condition_list:
        add_diagnostic_survey(cb, coverage=coverage, repetitions=repetitions, tsteps_btwn=interval,
                              target=target_group, start_day=start_days[0],
                              diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold,
                              nodes=nodes, positive_diagnosis_configs=msat_cfg,
                              IP_restrictions=ind_property_restrictions, NP_restrictions=node_property_restrictions,
                              pos_diag_IP_restrictions=IP_restrictions, trigger_condition_list=trigger_condition_list,
                              listening_duration=listening_duration, triggered_campaign_delay=triggered_campaign_delay)
    else:
        # MSAT controlled by MalariaDiagnostic campaign event rather than New_Diagnostic_Sensitivity
        for start_day in start_days:
            add_diagnostic_survey(cb, coverage=coverage, repetitions=repetitions, tsteps_btwn=interval,
                                  target=target_group, start_day=start_day,
                                  diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold,
                                  nodes=nodes, positive_diagnosis_configs=msat_cfg,
                                  IP_restrictions=ind_property_restrictions, NP_restrictions=node_property_restrictions,
                                  pos_diag_IP_restrictions=IP_restrictions)


def add_fMDA(cb, start_days, trigger_coverage, coverage, drug_configs, receiving_drugs_event, repetitions, interval,
             treatment_delay, diagnostic_threshold, fmda_radius, node_selection_type, nodes, expire_recent_drugs,
             node_property_restrictions, ind_property_restrictions, target_group, trigger_condition_list,
             listening_duration, triggered_campaign_delay):
    fmda_trigger = "Give_Drugs_fMDA"
    fmda_setup = [fmda_cfg(fmda_radius, node_selection_type, event_trigger=fmda_trigger)]

    if treatment_delay > 0:
        fmda_setup = [{"class": "DelayedIntervention",
                       "Delay_Distribution": "FIXED_DURATION",
                       "Delay_Period": treatment_delay,
                       "Actual_IndividualIntervention_Configs": fmda_setup
                       }]

    event_config = drug_configs + [receiving_drugs_event]
    if expire_recent_drugs:
        event_config.append(expire_recent_drugs)

    if trigger_condition_list:
        add_diagnostic_survey(cb, coverage=trigger_coverage, repetitions=1, tsteps_btwn=interval,
                              target=target_group, start_day=start_days[0],
                              diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold,
                              nodes=nodes, positive_diagnosis_configs=fmda_setup,
                              IP_restrictions=ind_property_restrictions, NP_restrictions=node_property_restrictions,
                              trigger_condition_list=trigger_condition_list,
                              listening_duration=listening_duration, triggered_campaign_delay=triggered_campaign_delay)
        fmda_distribute_drugs = {"Event_Name": "Distribute fMDA",
                                 "class": "CampaignEvent",
                                 "Start_Day": start_days[0],
                                 "Event_Coordinator_Config":
                                     {
                                         "class": "StandardInterventionDistributionEventCoordinator",
                                         "Intervention_Config": {
                                             "class": "NodeLevelHealthTriggeredIV",
                                             "Demographic_Coverage": coverage,
                                             "Duration": -1,# runs forever listening for fmda triggers of which there might be many.
                                             "Trigger_Condition_List": [fmda_trigger],
                                             "Actual_IndividualIntervention_Config": {
                                                 "Intervention_List": event_config,
                                                 "class": "MultiInterventionDistributor"
                                             }
                                         }
                                     },
                                 "Nodeset_Config": nodes
                                 }
        if node_property_restrictions:
            fmda_distribute_drugs['Event_Coordinator_Config']['Node_Property_Restrictions'] = node_property_restrictions

        if expire_recent_drugs:
            fmda_distribute_drugs['Event_Coordinator_Config']["Intervention_Config"][
                "Property_Restrictions_Within_Node"] = [{"DrugStatus": "None"}]
        cb.add_event(fmda_distribute_drugs)

    else:
        for start_day in start_days:
            # separate event for each repetition, otherwise RCD and fMDA can get entangled.
            for rep in range(repetitions):
                add_diagnostic_survey(cb, coverage=trigger_coverage, repetitions=1, tsteps_btwn=interval,
                                      target=target_group, start_day=start_day + interval * rep,
                                      diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold,
                                      nodes=nodes, positive_diagnosis_configs=fmda_setup,
                                      IP_restrictions=ind_property_restrictions,
                                      NP_restrictions=node_property_restrictions)
                fmda_distribute_drugs = {"Event_Name": "Distribute fMDA",
                                         "class": "CampaignEvent",
                                         "Start_Day": start_day + interval * rep + treatment_delay,
                                         "Event_Coordinator_Config":
                                             {
                                                 "class": "StandardInterventionDistributionEventCoordinator",
                                                 "Intervention_Config": {
                                                     "class": "NodeLevelHealthTriggeredIV",
                                                     "Demographic_Coverage": coverage,
                                                     "Duration": 2, #no confusion, but might as well just leave this 2 days.
                                                     "Trigger_Condition_List": [fmda_trigger],
                                                     "Actual_IndividualIntervention_Config": {
                                                         "Intervention_List": event_config,
                                                         "class": "MultiInterventionDistributor"
                                                     }
                                                 }
                                             },
                                         "Nodeset_Config": nodes
                                         }

                if node_property_restrictions:
                    fmda_distribute_drugs['Event_Coordinator_Config'][
                        'Node_Property_Restrictions'] = node_property_restrictions

                if expire_recent_drugs:
                    fmda_distribute_drugs['Event_Coordinator_Config']["Intervention_Config"][
                        "Property_Restrictions_Within_Node"] = [{"DrugStatus": "None"}]
                cb.add_event(fmda_distribute_drugs)


def add_rfMSAT(cb, start_day, coverage, drug_configs, receiving_drugs_event, interval, treatment_delay,
               trigger_coverage, diagnostic_threshold, fmda_radius, node_selection_type, snowballs, nodes,
               expire_recent_drugs, node_property_restrictions, ind_property_restrictions, trigger_condition_list,
                listening_duration, triggered_campaign_delay):
    fmda_setup = fmda_cfg(fmda_radius, node_selection_type) # no trigger used
    snowball_setup = [deepcopy(fmda_setup) for x in range(snowballs + 1)]
    snowball_trigger = 'Diagnostic_Survey_'
    snowball_setup[0]['Event_Trigger'] = snowball_trigger + "0"

    rcd_event = {"Event_Name": "Trigger RCD MSAT",
                 "class": "CampaignEvent",
                 "Start_Day": start_day,
                 "Event_Coordinator_Config":
                     {
                         "class": "StandardInterventionDistributionEventCoordinator",
                         "Intervention_Config": {
                             "class": "NodeLevelHealthTriggeredIV",
                             "Blackout_On_First_Occurrence": 1,
                             "Demographic_Coverage": trigger_coverage,
                             "Trigger_Condition_List": ["Received_Treatment"] + trigger_condition_list,
                             # triggered by successful health-seeking or a trigger that you give it.
                             "Duration": interval,  # interval argument indicates how long RCD will be implemented
                             "Actual_IndividualIntervention_Config": {
                                 "class": "DelayedIntervention",
                                 "Delay_Distribution": "FIXED_DURATION",
                                 "Delay_Period": 0,
                                 "Actual_IndividualIntervention_Configs": [snowball_setup[0]]
                             }
                         }
                     },
                 "Nodeset_Config": nodes}

    if node_property_restrictions:
        rcd_event['Event_Coordinator_Config']['Node_Property_Restrictions'] = node_property_restrictions

    cb.add_event(rcd_event)

    event_config = drug_configs + [receiving_drugs_event]
    IP_restrictions = []
    if expire_recent_drugs:
        event_config.append(expire_recent_drugs)
        IP_restrictions = [{"DrugStatus": "None"}]

    for snowball in range(snowballs):
        snowball_setup[snowball+1]['Event_Trigger'] = snowball_trigger + str(snowball+1)
        event_config = [snowball_setup[snowball+1], receiving_drugs_event] + drug_configs
        curr_trigger = snowball_trigger + str(snowball)
        add_diagnostic_survey(cb, coverage=coverage, start_day=start_day,
                              diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold, nodes=nodes,
                              trigger_condition_list=[curr_trigger], event_name='Snowball level ' + str(snowball),
                              positive_diagnosis_configs=event_config,
                              IP_restrictions=ind_property_restrictions, NP_restrictions=node_property_restrictions,
                              pos_diag_IP_restrictions=IP_restrictions)


def add_rfMDA(cb, start_day, coverage, drug_configs, receiving_drugs_event, interval, treatment_delay,
              trigger_coverage, fmda_radius, node_selection_type, nodes, expire_recent_drugs,
              node_property_restrictions, ind_property_restrictions, trigger_condition_list,
                listening_duration, triggered_campaign_delay=0):
    rfmda_trigger = "Give_Drugs_rfMDA"
    fmda_setup = fmda_cfg(fmda_radius, node_selection_type, event_trigger=rfmda_trigger)
    event_config = drug_configs + [receiving_drugs_event]
    if expire_recent_drugs:
        event_config.append(expire_recent_drugs)

    # distributes drugs to individuals broadcasting "Give_Drugs_rfMDA"
    # who is broadcasting is determined by other events
    # if campaign drugs change (less effective, different cocktail), then this event should have an expiration date.
    fmda_distribute_drugs = {"Event_Name": "Distribute fMDA",
                             "class": "CampaignEvent",
                             "Start_Day": start_day,
                             "Event_Coordinator_Config":
                                 {
                                     "class": "StandardInterventionDistributionEventCoordinator",
                                     "Intervention_Config": {
                                         "class": "NodeLevelHealthTriggeredIV",
                                         "Demographic_Coverage": coverage,
                                         "Duration": interval,
                                     # interval argument indicates how long RCD will be implemented
                                         "Trigger_Condition_List": [rfmda_trigger],
                                         "Actual_IndividualIntervention_Config": {
                                             "Intervention_List": event_config,
                                             "class": "MultiInterventionDistributor"
                                         }
                                     }
                                 },
                             "Nodeset_Config": nodes
                             }
    if node_property_restrictions:
        fmda_distribute_drugs['Event_Coordinator_Config']['Node_Property_Restrictions'] = node_property_restrictions
    if expire_recent_drugs:
        fmda_distribute_drugs['Event_Coordinator_Config']["Intervention_Config"][
            "Property_Restrictions_Within_Node"] = [{"DrugStatus": "None"}]

    rcd_event = {"Event_Name": "Trigger RCD MDA",
                 "class": "CampaignEvent",
                 "Start_Day": start_day,
                 "Event_Coordinator_Config":
                     {
                         "class": "StandardInterventionDistributionEventCoordinator",
                         "Intervention_Config": {
                             "class": "NodeLevelHealthTriggeredIV",
                             "Demographic_Coverage": trigger_coverage,
                             "Trigger_Condition_List": ["Received_Treatment"] + trigger_condition_list,  # triggered by successful health-seeking
                             "Duration": interval,  # interval argument indicates how long RCD will be implemented
                             "Actual_IndividualIntervention_Config": {
                                 "class": "DelayedIntervention",
                                 "Delay_Distribution": "FIXED_DURATION",
                                 "Delay_Period": triggered_campaign_delay,
                                 "Actual_IndividualIntervention_Configs": [fmda_setup]
                             }
                         }
                     },
                 "Nodeset_Config": nodes}

    if ind_property_restrictions:
        rcd_event['Event_Coordinator_Config']["Property_Restrictions_Within_Node"] = ind_property_restrictions
    if node_property_restrictions:
        rcd_event['Event_Coordinator_Config']['Node_Property_Restrictions'] = node_property_restrictions

    cb.add_event(rcd_event)
    cb.add_event(fmda_distribute_drugs)


# not implemented, do not update.
def add_borderscreen(cb, start_day, coverage, drug_configs, receiving_drugs_event, interval, treatment_delay,
                     trigger_coverage, diagnostic_threshold, fmda_radius, node_selection_type, snowballs, nodes,
                     expire_recent_drugs):
    # NOT TESTED; RECENT DRUG INELIGIBILITY ALSO NOT TESTED OR FULLY IMPLEMENTED

    fmda_setup = fmda_cfg(fmda_radius, node_selection_type)
    snowball_setup = [deepcopy(fmda_setup) for x in range(snowballs + 1)]

    IP_restrictions = []
    if expire_recent_drugs:
        IP_restrictions = [{"DrugStatus": "None"}]

    rcd_event = {"Event_Name": "Trigger Border Screen",
                 "class": "CampaignEvent",
                 "Start_Day": start_day,
                 "Event_Coordinator_Config":
                     {
                         "class": "StandardInterventionDistributionEventCoordinator",
                         "Intervention_Config": {
                             "class": "NodeLevelHealthTriggeredIV",
                             "Demographic_Coverage": trigger_coverage,
                             "Trigger_Condition_List": ["Emigrating"],  # triggered by successful health-seeking
                             "Duration": interval,  # interval argument indicates how long RCD will be implemented
                             "Actual_IndividualIntervention_Config": {
                                 "class": "DelayedIntervention",
                                 "Delay_Distribution": "FIXED_DURATION",
                                 "Delay_Period": treatment_delay,
                                 "Actual_IndividualIntervention_Configs": [snowball_setup[0]]
                             }
                         }
                     },
                 "Nodeset_Config": nodes}
    cb.add_event(rcd_event)
    add_diagnostic_survey(cb, coverage=coverage, start_day=start_day,
                          diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold, nodes=nodes,
                          trigger_condition_list=[snowball_setup[0]['Event_Trigger']],
                          event_name='Reactive MSAT level 0',
                          positive_diagnosis_configs=drug_configs + [receiving_drugs_event],
                          IP_restrictions=IP_restrictions)

    if snowballs > 0:
        for snowball in range(snowballs + 1):
            snowball_setup[snowball]['Event_Trigger'] = 'Diagnostic_Survey_' + str(snowball)
        for snowball in range(snowballs):
            curr_trigger = snowball_setup[snowball]['Event_Trigger']
            add_diagnostic_survey(cb, coverage=coverage, start_day=start_day,
                                  diagnostic_type='Other', diagnostic_threshold=diagnostic_threshold, nodes=nodes,
                                  trigger_condition_list=[curr_trigger], event_name='Snowball level ' + str(snowball),
                                  positive_diagnosis_configs=[snowball_setup[snowball + 1],
                                                              receiving_drugs_event] + drug_configs,
                                  IP_restrictions=IP_restrictions)


def fmda_cfg(fmda_type, node_selection_type='DISTANCE_ONLY', event_trigger='Give_Drugs'):
    # defaults to household-only unless fmda_type is a distance in km
    # options on Node_Selection_Type :
    # DISTANCE_ONLY: It will send the event to nodes that are within a given distance.
    # MIGRATION_NODES_ONLY: It will only send the event to nodes that the individual can migrate to.
    # DISTANCE_AND_MIGRATION: It will only send the even to migratable nodes that are within a given distance.
    # Migrateable nodes = Local and Regional.

    fmda = {
        "class": "BroadcastEventToOtherNodes",
        "Event_Trigger": event_trigger,
        "Include_My_Node": 1,
        "Node_Selection_Type": node_selection_type,
        "Max_Distance_To_Other_Nodes_Km": 0
    }

    if isinstance(fmda_type, int) or isinstance(fmda_type, float):
        fmda["Max_Distance_To_Other_Nodes_Km"] = fmda_type
    else:
        try:
            fmda["Max_Distance_To_Other_Nodes_Km"] = float(fmda_type)
        except ValueError:
            pass
    return fmda
